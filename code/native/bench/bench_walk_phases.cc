// B1 -- phase profile of the campaign's inner loop, by ablation.
//
// The walk is written once and instantiated at four levels; each level adds one
// phase, so consecutive differences are the phases:
//
//   L0  skeleton   the DFS itself: prime lookups, bound tests, divisor counts
//   L1  + push     Spine::push / pop -- multiply, canonize, sort, unique
//   L2  + nu       the nu(z) = Re(z^2) extraction loop at emission
//   L3  + meet     meet<trunc_key_t> with the D-dispatch, i.e. the campaign
//
// Ablation rather than sampling, for a reason worth keeping: the tree's shape
// depends only on primes, bounds and exponent counts -- never on points -- so
// removing `push` traverses *exactly* the same tree. The bench asserts that by
// requiring all four levels to report the same norm count, which also makes it
// self-checking against the fixtures (11,083 / 131,437 / 1,464,114).
//
//   ./bench_walk_phases <bound> <repeats>
#include "bench/harness.h"
#include "cpu/core/precision.h"
#include "cpu/factorbase/factorbase.h"
#include "cpu/gaussian/gaussian.h"
#include "cpu/generate/generate.h"
#include "cpu/meet/meet.h"

#include <cstdio>
#include <cstdlib>
#include <span>
#include <vector>

using namespace pte::cpu;
using pte::bench::Result;
using pte::bench::Suite;

namespace {

enum Level { kSkeleton = 0, kPush = 1, kNu = 2, kMeet = 3 };

struct Counters {
  std::size_t norms = 0;
  std::size_t points = 0;
  std::size_t partitions = 0; // sum of C(D,4)*3, the direct-meet work metric
  std::size_t hits = 0;
  std::uint64_t checksum = 0;
};

std::size_t QuadPartitions(std::size_t d) {
  return d < 4 ? 0 : d * (d - 1) * (d - 2) * (d - 3) / 8;
}

// A transcription of Walker::walk. Deliberately a copy rather than a template
// parameter on the shipped Walker: the production header should not carry a
// knob that exists only for measurement.
template <Level L> class PhaseWalker {
public:
  explicit PhaseWalker(const Generator &g)
      : g_(&g), spine_(g.base()), nu_(kMaxPointCount) {}

  void walk(std::size_t segment, Counters &c) {
    const norm_t bound = g_->bound();
    const FactorBase &base = g_->base();
    const auto top =
        static_cast<norm_t>(base.prime(static_cast<factor_index_t>(segment)));
    if (top > bound)
      return;

    const auto emit = [&](norm_t v) {
      ++c.norms;
      if constexpr (L == kSkeleton) {
        c.checksum += static_cast<std::uint64_t>(v);
      } else {
        const std::span<const Gauss> points = spine_.points();
        c.points += points.size();
        c.partitions += QuadPartitions(points.size());
        if constexpr (L == kPush) {
          c.checksum += static_cast<std::uint64_t>(points[0].re);
        } else {
          for (std::size_t k = 0; k < points.size(); ++k)
            nu_[k] = static_cast<norm_t>(points[k].nu());
          if constexpr (L == kNu) {
            c.checksum += static_cast<std::uint64_t>(nu_[0]);
          } else {
            meet<trunc_key_t>(static_cast<int>(points.size()), nu_.data(),
                              scratch_, [&](const E4Sol &) { ++c.hits; });
          }
        }
      }
    };

    const auto push = [&](std::size_t i) {
      if constexpr (L != kSkeleton)
        spine_.push(static_cast<factor_index_t>(i));
    };
    const auto pop = [&] {
      if constexpr (L != kSkeleton)
        spine_.pop();
    };

    const auto dfs = [&](auto &&self, std::size_t start, norm_t value,
                         std::size_t divisors) -> void {
      for (std::size_t i = start; i < segment; ++i) {
        const auto p =
            static_cast<norm_t>(base.prime(static_cast<factor_index_t>(i)));
        if (value > bound / p)
          break;
        norm_t v = value * p;
        std::size_t exponent = 1;
        while (true) {
          push(i);
          const std::size_t d = divisors * (exponent + 1);
          if (d >= 8)
            emit(v);
          self(self, i + 1, v, d);
          if (v > bound / p)
            break;
          v *= p;
          ++exponent;
        }
        for (std::size_t e = 0; e < exponent; ++e)
          pop();
      }
    };

    norm_t v = top;
    std::size_t exponent = 1;
    while (true) {
      push(segment);
      const std::size_t d = exponent + 1;
      if (d >= 8)
        emit(v);
      dfs(dfs, 0, v, d);
      if (v > bound / top)
        break;
      v *= top;
      ++exponent;
    }
    for (std::size_t e = 0; e < exponent; ++e)
      pop();
  }

private:
  const Generator *g_;
  Spine spine_;
  std::vector<norm_t> nu_;
  MeetScratch<trunc_key_t> scratch_;
};

template <Level L> void RunLevel(const Generator &g, Counters &c) {
  c = Counters{};
  PhaseWalker<L> w(g);
  for (std::size_t s = 0; s < g.segment_count(); ++s)
    w.walk(s, c);
}

} // namespace

int main(int argc, char **argv) {
  const norm_t bound = argc > 1 ? std::atoll(argv[1]) : 100'000'000;
  const int repeats = argc > 2 ? std::atoi(argv[2]) : 5;

  std::fprintf(stderr, "building factor base for B = %lld ...\n",
               static_cast<long long>(bound));
  const Generator g(bound);

  Counters c0, c1, c2, c3;
  RunLevel<kMeet>(g, c3); // once up front, for the unit count and the checks

  Suite suite("B1 walk phases (cumulative levels)",
              static_cast<double>(c3.norms), "norm");
  suite.Note("  " + std::to_string(c3.points) + " points, " +
             std::to_string(c3.partitions) + " quad-partitions, " +
             std::to_string(c3.hits) + " meet hits");
  suite.Add("L0 skeleton", [&] { RunLevel<kSkeleton>(g, c0); });
  suite.Add("L1 + push", [&] { RunLevel<kPush>(g, c1); });
  suite.Add("L2 + nu", [&] { RunLevel<kNu>(g, c2); });
  suite.Add("L3 + meet (campaign)", [&] { RunLevel<kMeet>(g, c3); });
  const std::vector<Result> &r = suite.Run(repeats);

  // Self-check: the ablation is only valid if every level walked the same tree.
  if (c0.norms != c3.norms || c1.norms != c3.norms || c2.norms != c3.norms) {
    std::fprintf(stderr,
                 "ABLATION BROKEN: norm counts differ (%zu %zu %zu %zu)\n",
                 c0.norms, c1.norms, c2.norms, c3.norms);
    return 1;
  }

  const double n = static_cast<double>(c3.norms);
  const double ns = 1e9 / n;
  const bool have = r[0].counters.available;
  std::printf("\n  %-24s %9s %9s %8s", "phase (difference)", "ns/norm",
              "cyc/norm", "share");
  if (have)
    std::printf(" %9s", "brmiss/n");
  std::printf("\n");

  const char *names[] = {"DFS skeleton", "Spine::push", "nu extraction", "meet"};
  const double total = r[3].min_s;
  for (int i = 0; i < 4; ++i) {
    const double dt = i == 0 ? r[0].min_s : r[i].min_s - r[i - 1].min_s;
    const double dc =
        i == 0 ? static_cast<double>(r[0].counters.cycles)
               : static_cast<double>(r[i].counters.cycles) -
                     static_cast<double>(r[i - 1].counters.cycles);
    const double db =
        i == 0 ? static_cast<double>(r[0].counters.branch_misses)
               : static_cast<double>(r[i].counters.branch_misses) -
                     static_cast<double>(r[i - 1].counters.branch_misses);
    std::printf("  %-24s %9.1f %9.1f %7.1f%%", names[i], dt * ns,
                have ? dc / n : 0.0, 100.0 * dt / total);
    if (have)
      std::printf(" %9.4f", db / n);
    std::printf("\n");
  }

  if (c3.partitions > 0) {
    const double dt = r[3].min_s - r[2].min_s;
    std::printf("\n  meet: %.2f ns / quad-partition over %zu partitions "
                "(%.1f per norm)\n",
                dt * 1e9 / static_cast<double>(c3.partitions), c3.partitions,
                static_cast<double>(c3.partitions) / n);
  }
  std::fprintf(stderr, "sink %llu\n",
               (unsigned long long)(c0.checksum + c1.checksum + c2.checksum +
                                    c3.checksum + suite.sink()));
  return 0;
}
