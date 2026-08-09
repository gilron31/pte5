#include "cpu/campaign/campaign.h"

#include "cpu/confirm/confirm.h"
#include "cpu/enumerate/enumerate.h"
#include "cpu/meet/meet.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unistd.h>
#include <utility>
#include <vector>

namespace pte::cpu {
namespace {

using Clock = std::chrono::steady_clock;

double SecondsSince(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

// Adaptive batching bounds. The opening batch is small because the opening
// segments are the expensive ones; the cap keeps the checkpoint from drifting
// arbitrarily far out on a fast machine.
constexpr std::size_t kInitialBatch = 64;
constexpr std::size_t kMaxBatch = 1u << 20;

// The checkpoint. Deliberately a plain key=value text file: it is the one
// artifact a human reads when a multi-hour run has to be restarted, and a
// binary format would make "how far did it get" a debugging exercise.
struct Checkpoint {
  int version = 1;
  norm_t bound = 0;
  std::size_t segments_done = 0;
  std::size_t catalog_bytes = 0;
  CampaignStats stats;

  std::string Serialize() const {
    std::ostringstream os;
    os << "version=" << version << "\n"
       << "bound=" << bound << "\n"
       << "segments_done=" << segments_done << "\n"
       << "catalog_bytes=" << catalog_bytes << "\n"
       << "norms=" << stats.norms << "\n"
       << "candidates=" << stats.candidates << "\n"
       << "filtered=" << stats.filtered << "\n"
       << "imprimitive=" << stats.imprimitive << "\n"
       << "uncertified=" << stats.uncertified << "\n"
       << "solutions=" << stats.solutions << "\n";
    return os.str();
  }

  static Checkpoint Parse(const std::string &path) {
    std::ifstream in(path);
    if (!in)
      throw std::runtime_error("campaign: cannot read checkpoint " + path);
    Checkpoint c;
    std::string line;
    while (std::getline(in, line)) {
      const auto eq = line.find('=');
      if (eq == std::string::npos)
        continue;
      const std::string key = line.substr(0, eq);
      const std::string value = line.substr(eq + 1);
      if (key == "version")
        c.version = std::stoi(value);
      else if (key == "bound")
        c.bound = std::stoll(value);
      else if (key == "segments_done")
        c.segments_done = std::stoull(value);
      else if (key == "catalog_bytes")
        c.catalog_bytes = std::stoull(value);
      else if (key == "norms")
        c.stats.norms = std::stoull(value);
      else if (key == "candidates")
        c.stats.candidates = std::stoull(value);
      else if (key == "filtered")
        c.stats.filtered = std::stoull(value);
      else if (key == "imprimitive")
        c.stats.imprimitive = std::stoull(value);
      else if (key == "uncertified")
        c.stats.uncertified = std::stoull(value);
      else if (key == "solutions")
        c.stats.solutions = std::stoull(value);
    }
    if (c.version != 1)
      throw std::runtime_error("campaign: checkpoint " + path + " is version " +
                               std::to_string(c.version) + ", expected 1");
    return c;
  }
};

// Replace `path` atomically. The rename is the commit point; without the fsync
// before it, a power loss can leave the new name pointing at empty content.
void WriteAtomically(const std::string &path, const std::string &content) {
  const std::string tmp = path + ".tmp";
  {
    std::FILE *f = std::fopen(tmp.c_str(), "wb");
    if (f == nullptr)
      throw std::runtime_error("campaign: cannot write " + tmp);
    const bool ok =
        std::fwrite(content.data(), 1, content.size(), f) == content.size();
    std::fflush(f);
    ::fsync(::fileno(f));
    std::fclose(f);
    if (!ok)
      throw std::runtime_error("campaign: short write to " + tmp);
  }
  std::filesystem::rename(tmp, path);
}

// Per-thread accumulators, merged after the batch joins.
struct ThreadState {
  Walker walker;
  MeetScratch<trunc_key_t> scratch;
  std::vector<factor_t> primes; // reused across candidates
  CampaignStats stats;

  explicit ThreadState(const Generator &g) : walker(g) {}
};

// One catalog row. nu is in canonical solution order, so the row is a function
// of the solution rather than of the order the meet happened to emit it in.
std::string FormatRow(norm_t l1, const std::span<const norm_t> nu,
                      const E4Sol &canonical, const Certificate &cert) {
  std::ostringstream os;
  os << l1 << ',' << nu[static_cast<std::size_t>(canonical.a)] << ','
     << nu[static_cast<std::size_t>(canonical.b)] << ','
     << nu[static_cast<std::size_t>(canonical.c)] << ','
     << nu[static_cast<std::size_t>(canonical.d)] << ','
     << cert.tower.l2.get_str() << ',' << cert.tower.l3.get_str() << ','
     << cert.tower.sqrt_l4.get_str() << ',' << (cert.valid() ? 1 : 0) << '\n';
  return os.str();
}

} // namespace

const char *Campaign::header() {
  return "l1,nu0,nu1,nu2,nu3,l2,l3,sqrt_l4,certified";
}

namespace {

// Builds the Generator, loading the factor base from an artifact when one is
// given. Free function because it has to run inside Campaign's member-init
// list.
Generator MakeGenerator(const CampaignOptions &opts, double &seconds) {
  const auto start = Clock::now();
  Generator g = opts.factor_base.empty()
                    ? Generator(opts.bound)
                    : Generator(opts.bound,
                                FactorBase::Load(opts.factor_base,
                                                 static_cast<factor_t>(
                                                     opts.bound /
                                                     kSmallestCofactor)));
  seconds = SecondsSince(start);
  return g;
}

} // namespace

Campaign::Campaign(CampaignOptions opts)
    : opts_(std::move(opts)),
      generator_(MakeGenerator(opts_, base_seconds_)) {
  if (opts_.checkpoint.empty())
    opts_.checkpoint = opts_.out + ".ckpt";
  if (opts_.checkpoint_seconds <= 0)
    throw std::invalid_argument(
        "campaign: checkpoint interval must be positive");
  if (opts_.threads <= 0) {
    const unsigned hw = std::thread::hardware_concurrency();
    opts_.threads = static_cast<int>(hw == 0 ? 1u : hw);
  }
}

CampaignStats Campaign::run() {
  const std::size_t segments = generator_.segment_count();
  CampaignStats stats;
  stats.base_seconds = base_seconds_;

  std::size_t start_segment = 0;
  std::size_t catalog_bytes = 0;

  if (opts_.resume && std::filesystem::exists(opts_.checkpoint)) {
    const Checkpoint c = Checkpoint::Parse(opts_.checkpoint);
    if (c.bound != opts_.bound) {
      throw std::runtime_error(
          "campaign: checkpoint is for bound " + std::to_string(c.bound) +
          ", not " + std::to_string(opts_.bound) +
          "; resuming would merge two different searches into one catalog");
    }
    start_segment = c.segments_done;
    catalog_bytes = c.catalog_bytes;
    stats = c.stats;
    // segments_done is a top-level checkpoint key, not part of the stats block,
    // so it has to be restored explicitly. Without this a resume of an
    // *already finished* campaign reports INCOMPLETE -- the batch loop never
    // runs, so nothing ever sets it -- and an operator re-runs forever.
    stats.segments_done = c.segments_done;
    stats.base_seconds = base_seconds_;

    // Discard anything written after the last durable checkpoint. Without this,
    // a batch that reached the catalog but not the checkpoint would be counted
    // twice.
    if (std::filesystem::exists(opts_.out))
      std::filesystem::resize_file(opts_.out, catalog_bytes);
  }

  const bool fresh = start_segment == 0;
  std::FILE *catalog = std::fopen(opts_.out.c_str(), fresh ? "wb" : "r+b");
  if (catalog == nullptr)
    throw std::runtime_error("campaign: cannot open catalog " + opts_.out);
  if (fresh) {
    const std::string head = std::string(header()) + "\n";
    std::fwrite(head.data(), 1, head.size(), catalog);
    catalog_bytes = head.size();
  } else if (std::fseek(catalog, static_cast<long>(catalog_bytes), SEEK_SET) !=
             0) {
    std::fclose(catalog);
    throw std::runtime_error("campaign: cannot seek the catalog to its "
                             "checkpointed length");
  }

  const bool progress = !opts_.quiet && ::isatty(::fileno(stdout)) != 0;
  const auto walk_start = Clock::now();
  std::size_t batches = 0;

  // One Walker per thread for the whole run, not per batch: the Spine's buffer
  // is ~68 KB, and a B = 10^10 campaign has a thousand batches.
  std::vector<ThreadState> states;
  states.reserve(static_cast<std::size_t>(opts_.threads));
  for (int t = 0; t < opts_.threads; ++t)
    states.emplace_back(generator_);

  // Adaptive batching starts deliberately small: the first segments are the
  // expensive ones, so a large opening batch is exactly the wrong guess.
  const bool adaptive = opts_.batch == 0;
  std::size_t batch = adaptive ? kInitialBatch : opts_.batch;

  for (std::size_t begin = start_segment; begin < segments;) {
    if (opts_.stop_after_batches != 0 && batches == opts_.stop_after_batches)
      break;
    const std::size_t end = std::min(begin + batch, segments);
    const auto batch_start = Clock::now();

    // Fixed before the threads start, so concurrent writes land in distinct
    // elements and nothing reallocates.
    std::vector<std::string> rows(end - begin);
    std::atomic<std::size_t> next{begin};

    const auto work = [&](ThreadState &ts) {
      for (std::size_t s = next.fetch_add(1); s < end;
           s = next.fetch_add(1)) {
        std::string &sink = rows[s - begin];
        ts.walker.walk(s, [&](const NormView &norm) {
          ++ts.stats.norms;
          meet<trunc_key_t>(
              norm.D(), norm.nu.data(), ts.scratch, [&](const E4Sol &hit) {
                ++ts.stats.candidates;
                // Layer 2: the truncated key is a filter above kExactKeyBound,
                // so this is what turns a candidate into a result.
                if (!confirm_candidate(norm.nu, hit)) {
                  ++ts.stats.filtered;
                  return;
                }
                ts.primes.clear();
                for (const factor_index_t i : norm.factors)
                  ts.primes.push_back(generator_.base().prime(i));
                if (!is_primitive(norm.points, ts.primes, hit)) {
                  ++ts.stats.imprimitive;
                  return;
                }
                const Certificate cert =
                    certify(quadruple(norm.points, hit));
                if (!cert.valid())
                  ++ts.stats.uncertified;
                ++ts.stats.solutions;
                sink += FormatRow(norm.l1, norm.nu, hit.canonical(), cert);
              });
        });
      }
    };

    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(opts_.threads) - 1);
    for (int t = 1; t < opts_.threads; ++t)
      pool.emplace_back(work, std::ref(states[static_cast<std::size_t>(t)]));
    work(states[0]); // this thread pulls too rather than idling
    for (std::thread &th : pool)
      th.join();

    // Drain rather than accumulate, since the states outlive the batch.
    for (ThreadState &ts : states) {
      stats.norms += ts.stats.norms;
      stats.candidates += ts.stats.candidates;
      stats.filtered += ts.stats.filtered;
      stats.imprimitive += ts.stats.imprimitive;
      stats.uncertified += ts.stats.uncertified;
      stats.solutions += ts.stats.solutions;
      ts.stats = CampaignStats{};
    }

    // In segment order, so the catalog does not depend on the thread count.
    for (const std::string &row : rows) {
      if (row.empty())
        continue;
      if (std::fwrite(row.data(), 1, row.size(), catalog) != row.size()) {
        std::fclose(catalog);
        throw std::runtime_error("campaign: short write to " + opts_.out);
      }
      catalog_bytes += row.size();
    }
    std::fflush(catalog);
    ::fsync(::fileno(catalog));

    stats.segments_done = end;
    Checkpoint ck;
    ck.bound = opts_.bound;
    ck.segments_done = end;
    ck.catalog_bytes = catalog_bytes;
    ck.stats = stats;
    WriteAtomically(opts_.checkpoint, ck.Serialize());

    ++batches;
    begin = end;

    // Rescale towards the target interval. Clamped per step so one anomalous
    // batch -- a page-cache stall, another process -- cannot swing the size by
    // orders of magnitude.
    if (adaptive) {
      const double took = SecondsSince(batch_start);
      double scale = took > 0 ? opts_.checkpoint_seconds / took : 4.0;
      scale = std::clamp(scale, 0.25, 4.0);
      const auto scaled =
          static_cast<std::size_t>(static_cast<double>(batch) * scale);
      batch = std::clamp<std::size_t>(scaled, 1, kMaxBatch);
    }

    // The in-place progress line is for a human watching a terminal. Redirected
    // to a file or a pipe the carriage returns do not overwrite anything, and a
    // thousand-batch run turns into tens of kilobytes of noise around the
    // summary -- so off unless stdout is a terminal.
    if (progress) {
      std::printf("\r  segments %zu/%zu  norms %zu  solutions %zu   ", end,
                  segments, stats.norms, stats.solutions);
      std::fflush(stdout);
    }
  }

  std::fclose(catalog);
  stats.walk_seconds = SecondsSince(walk_start);
  stats.complete = stats.segments_done >= segments;
  if (progress)
    std::printf("\n");
  return stats;
}

} // namespace pte::cpu
