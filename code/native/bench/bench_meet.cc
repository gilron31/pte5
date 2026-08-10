// B3 -- the meet, decomposed.
//
// Why this one matters most: bench_walk_phases shows the meet at 42% of the
// inner loop at B = 1e7, 50% at 1e8 and 57% at 1e9, growing while Spine::push
// stays flat. At the design bound it is the campaign. So this bench asks three
// questions the phase profile cannot:
//
//   1. Per D, where is the direct/hashed crossover, in cycles rather than
//      wall clock, and with branch-miss counts attached?
//   2. Inside meet_hashed -- which carries 74% of direct-meet work at 1e8 --
//      is the cost the key generation, the sort, or the scan?
//   3. How far is any of it from the ALU ceiling?
//
// The corpus is every admissible norm's nu array at the given bound, harvested
// once and grouped by D, so each measurement runs on the real distribution
// rather than a synthetic one. drafts/implementation_design.md 7.1 records what
// happens when you forget that: a fixed-D microbenchmark ranks hashed *above*
// meet_direct<8>, and the mixed corpus reverses it.
//
//   ./bench_meet <bound> <repeats>
#include "bench/ceiling.h"
#include "bench/harness.h"
#include "cpu/core/precision.h"
#include "cpu/core/solution.h"
#include "cpu/generate/generate.h"
#include "cpu/meet/meet.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <type_traits>
#include <vector>

using namespace pte::cpu;
using pte::bench::Ceiling;
using pte::bench::Machine;
using pte::bench::PerfCounters;
using pte::bench::Stopwatch;
using pte::bench::Suite;

namespace {

// ---------------------------------------------------------------------------
// Corpus
// ---------------------------------------------------------------------------

struct Corpus {
  std::vector<norm_t> nu;        // every norm's values, packed end to end
  std::vector<std::size_t> at;   // norm i owns nu[at[i] .. at[i] + d[i])
  std::vector<int> d;
  std::vector<std::vector<std::size_t>> by_d; // indices into `at`, per D
  std::size_t partitions = 0;    // sum of C(D,4)*3
  int max_d = 0;
};

std::size_t QuadPartitions(int d) {
  return d < 4 ? 0
               : static_cast<std::size_t>(d) * (d - 1) * (d - 2) * (d - 3) / 8;
}

// Cycles for one norm's direct meet, if nothing ever stalled.
//
// The **hoist** is D independent 64-bit squarings. Zen 4 retires one imul per
// cycle, so this is D cycles flat and it does *not* divide by the ALU width. At
// D = 4 that is 4 cycles against 1.5 for the three partitions: the modal norm
// spends most of its meet squaring, not comparing.
//
// The **table loop** is where a first draft of this model went wrong, and the
// error is worth keeping written down because the corrected version is the
// actionable result. It is tempting to call the loop ALU-bound at 1 add + 1 cmp
// per partition with sq[] in registers. But `partitions<D>()` is a constexpr
// *array*, not unrolled code, so each iteration:
//
//     loads one 4-byte Partition {a,b,c,d}                        1 load
//     indexes sq[] four times with those runtime indices          4 loads
//     one add, one compare, plus loop increment and branch      ~3 ALU ops
//
// Five loads against three load ports is 1.67 cycles; three ALU ops against
// four pipes is 0.75. **The direct meet is load-bound, and by more than 2x.**
// Which says something specific: the way past this ceiling is not wider adds,
// it is turning those indices from data into code -- compile-time unrolling, or
// SIMD shuffles over a fixed permutation -- so the four gathers disappear.
double DirectCeilingCycles(int d) {
  const double hoist = static_cast<double>(d) / Machine::kImul64PerCycle;
  const double p = static_cast<double>(QuadPartitions(d));
  const double load_bound = p * 5.0 / Machine::kLoadPorts;
  const double alu_bound = p * 3.0 / Machine::kIntAluPipes;
  return hoist + std::max(load_bound, alu_bound);
}

Corpus Harvest(const Generator &g) {
  Corpus c;
  g.for_each_norm([&](const NormView &n) {
    c.at.push_back(c.nu.size());
    c.d.push_back(n.D());
    c.max_d = std::max(c.max_d, n.D());
    c.partitions += QuadPartitions(n.D());
    for (const norm_t v : n.nu)
      c.nu.push_back(v);
  });
  c.by_d.resize(static_cast<std::size_t>(c.max_d) + 1);
  for (std::size_t i = 0; i < c.at.size(); ++i)
    c.by_d[static_cast<std::size_t>(c.d[i])].push_back(i);
  return c;
}

// ---------------------------------------------------------------------------
// meet_hashed, cut into stages
// ---------------------------------------------------------------------------
//
// A transcription of the shipped meet_hashed, stopped early. Transcription is a
// hazard -- it can drift from the original -- so kFull is run against the real
// meet_hashed and the hit counts are required to agree before any number below
// is printed.

enum class Stage { kHoist, kKeys, kSort, kFull };

template <Stage S, typename Sink>
void HashedUpTo(int D, const norm_t *nu, MeetScratch<trunc_key_t> &scratch,
                Sink &&sink) {
  if (D < 4)
    return;
  auto &sq = scratch.sq;
  auto &pairs = scratch.pairs;
  sq.resize(static_cast<std::size_t>(D));
  for (int i = 0; i < D; ++i) {
    const auto v = static_cast<trunc_key_t>(nu[i]);
    sq[static_cast<std::size_t>(i)] = v * v;
  }
  if constexpr (S == Stage::kHoist)
    return;

  pairs.clear();
  pairs.reserve(static_cast<std::size_t>(D) * (D - 1) / 2);
  for (int i = 0; i < D; ++i)
    for (int j = i + 1; j < D; ++j)
      pairs.push_back(
          {sq[static_cast<std::size_t>(i)] + sq[static_cast<std::size_t>(j)],
           static_cast<std::uint16_t>(i), static_cast<std::uint16_t>(j)});
  if constexpr (S == Stage::kKeys)
    return;

  std::sort(pairs.begin(), pairs.end(),
            [](const auto &x, const auto &y) { return x.key < y.key; });
  if constexpr (S == Stage::kSort)
    return;

  for (std::size_t i = 0; i < pairs.size();) {
    std::size_t j = i + 1;
    while (j < pairs.size() && pairs[j].key == pairs[i].key)
      ++j;
    for (std::size_t p = i; p + 1 < j; ++p)
      for (std::size_t q = p + 1; q < j; ++q) {
        const auto &u = pairs[p];
        const auto &v = pairs[q];
        if (u.a == v.a || u.a == v.b || u.b == v.a || u.b == v.b)
          continue;
        sink(E4Sol{u.a, u.b, v.a, v.b});
      }
    i = j;
  }
}

// ---------------------------------------------------------------------------
// Compile-time D dispatch for the per-D sweep
// ---------------------------------------------------------------------------

// meet_direct needs D at compile time, so the sweep instantiates a range and
// looks the right one up. 16 is where the table stops being cheap:
// partition_count(16) = 5,460 entries, against 14,535 at D = 20.
inline constexpr int kMaxSweepD = 16;

template <int D, typename Sink>
void DirectAt(int d, const norm_t *nu, Sink &&sink) {
  if constexpr (D > kMaxSweepD) {
    (void)d;
    (void)nu;
    (void)sink;
  } else {
    if (d == D)
      meet_direct<D, trunc_key_t>(nu, sink);
    else
      DirectAt<D + 1>(d, nu, sink);
  }
}

struct Timing {
  double seconds = 0;
  PerfCounters::Reading counters;
};

template <typename F> Timing TimeIt(int repeats, F &&f) {
  Timing best;
  best.seconds = 1e300;
  for (int r = 0; r < repeats; ++r) {
    PerfCounters pc;
    Stopwatch sw;
    pc.Start();
    sw.Start();
    f();
    const double s = sw.CpuSeconds();
    const PerfCounters::Reading cr = pc.Stop();
    if (s < best.seconds) {
      best.seconds = s;
      best.counters = cr;
    }
  }
  return best;
}

} // namespace

int main(int argc, char **argv) {
  const norm_t bound = argc > 1 ? std::atoll(argv[1]) : 10'000'000;
  const int repeats = argc > 2 ? std::atoi(argv[2]) : 7;

  std::fprintf(stderr, "harvesting corpus at B = %lld ...\n",
               static_cast<long long>(bound));
  const Generator g(bound);
  const Corpus c = Harvest(g);
  const auto norms = static_cast<double>(c.at.size());
  if (c.at.empty()) {
    std::printf("no admissible norms below %lld\n",
                static_cast<long long>(bound));
    return 0;
  }

  // -------------------------------------------------------------------------
  // Correctness gate: the transcription must agree with the shipped function.
  // -------------------------------------------------------------------------
  {
    std::size_t mine = 0, theirs = 0;
    MeetScratch<trunc_key_t> s1, s2;
    for (std::size_t i = 0; i < c.at.size(); ++i) {
      if (c.d[i] <= 8)
        continue; // meet_hashed is only reached above the tables
      HashedUpTo<Stage::kFull>(c.d[i], c.nu.data() + c.at[i], s1,
                               [&](const E4Sol &) { ++mine; });
      meet_hashed<trunc_key_t>(c.d[i], c.nu.data() + c.at[i], s2,
                               [&](const E4Sol &) { ++theirs; });
    }
    if (mine != theirs) {
      std::fprintf(stderr,
                   "TRANSCRIPTION DRIFT: staged meet_hashed found %zu hits, "
                   "the shipped one %zu -- numbers below would be lies\n",
                   mine, theirs);
      return 1;
    }
    std::fprintf(stderr, "staged meet_hashed agrees with the shipped one (%zu "
                         "hits on D > 8)\n",
                 theirs);
  }

  std::size_t hits = 0;
  MeetScratch<trunc_key_t> scratch;
  const auto sink = [&](const E4Sol &) { ++hits; };

  // -------------------------------------------------------------------------
  // 1. Dispatch, over the whole corpus
  // -------------------------------------------------------------------------
  {
    Suite s("B3.1 meet dispatch over every admissible norm", norms, "norm");
    s.Note("  " + std::to_string(c.partitions) + " quad-partitions (" +
           std::to_string(static_cast<long long>(
               static_cast<double>(c.partitions) / norms)) +
           "/norm), max D = " + std::to_string(c.max_d));
    double ceil_cycles = 0;
    for (const int d : c.d)
      ceil_cycles += DirectCeilingCycles(d);
    s.SetCeiling(Ceiling{ceil_cycles / norms,
                         "D squarings at 1 imul/cyc, plus the table loop at "
                         "5 loads / 3 ports per quad-partition (load-bound)"});
    s.Add("meet<> (shipped dispatch)", [&] {
      for (std::size_t i = 0; i < c.at.size(); ++i)
        meet<trunc_key_t>(c.d[i], c.nu.data() + c.at[i], scratch, sink);
    });
    s.Add("all hashed", [&] {
      for (std::size_t i = 0; i < c.at.size(); ++i)
        meet_hashed<trunc_key_t>(c.d[i], c.nu.data() + c.at[i], scratch, sink);
    });
    s.Add("all direct (D <= 16 only)", [&] {
      for (std::size_t i = 0; i < c.at.size(); ++i) {
        if (c.d[i] <= kMaxSweepD)
          DirectAt<4>(c.d[i], c.nu.data() + c.at[i], sink);
        else
          meet_hashed<trunc_key_t>(c.d[i], c.nu.data() + c.at[i], scratch,
                                   sink);
      }
    });
    s.Run(repeats);
  }

  // -------------------------------------------------------------------------
  // 2. meet_hashed, by stage, on the D > 8 tail only
  // -------------------------------------------------------------------------
  {
    std::vector<std::size_t> tail;
    std::size_t tail_keys = 0;
    for (std::size_t i = 0; i < c.at.size(); ++i)
      if (c.d[i] > 8) {
        tail.push_back(i);
        tail_keys += static_cast<std::size_t>(c.d[i]) * (c.d[i] - 1) / 2;
      }
    if (!tail.empty()) {
      const auto tn = static_cast<double>(tail.size());
      Suite s("B3.2 meet_hashed by stage, D > 8 tail (cumulative)", tn, "norm");
      s.Note("  " + std::to_string(tail.size()) + " norms (" +
             std::to_string(100.0 * tn / norms).substr(0, 4) +
             "% of corpus), " + std::to_string(tail_keys) + " pair keys (" +
             std::to_string(static_cast<long long>(
                 static_cast<double>(tail_keys) / tn)) +
             "/norm)");
      const auto run = [&](auto stage_tag) {
        constexpr Stage St = decltype(stage_tag)::value;
        for (const std::size_t i : tail)
          HashedUpTo<St>(c.d[i], c.nu.data() + c.at[i], scratch, sink);
      };
      s.Add("S0 hoist sq[] only", [&] {
        run(std::integral_constant<Stage, Stage::kHoist>{});
      });
      s.Add("S1 + generate D(D-1)/2 keys", [&] {
        run(std::integral_constant<Stage, Stage::kKeys>{});
      });
      s.Add("S2 + sort", [&] {
        run(std::integral_constant<Stage, Stage::kSort>{});
      });
      s.Add("S3 + scan and emit (full)", [&] {
        run(std::integral_constant<Stage, Stage::kFull>{});
      });
      const auto &r = s.Run(repeats);

      std::printf("\n  %-30s %9s %9s %8s\n", "stage (difference)", "ns/norm",
                  "cyc/norm", "share");
      const char *names[] = {"hoist sq[]", "key generation", "sort",
                             "scan + emit"};
      for (int i = 0; i < 4; ++i) {
        const double dt =
            i == 0 ? r[0].min_s : r[i].min_s - r[i - 1].min_s;
        const double dc =
            i == 0 ? static_cast<double>(r[0].counters.cycles)
                   : static_cast<double>(r[i].counters.cycles) -
                         static_cast<double>(r[i - 1].counters.cycles);
        std::printf("  %-30s %9.1f %9.1f %7.1f%%\n", names[i], dt * 1e9 / tn,
                    dc / tn, 100.0 * dt / r[3].min_s);
      }
      std::printf("  (sort is over %.0f keys/norm; a negative difference means "
                  "the stages are within noise of each other)\n",
                  static_cast<double>(tail_keys) / tn);
    }
  }

  // -------------------------------------------------------------------------
  // 3. Per-D crossover
  // -------------------------------------------------------------------------
  std::printf("\n=== B3.3 direct vs hashed, per D, on the real corpus ===\n");
  std::printf("  %3s %9s %10s %10s %8s %10s %10s %9s\n", "D", "norms",
              "direct ns", "hashed ns", "winner", "direct cyc", "hashed cyc",
              "%ceil dir");
  for (int d = 4; d <= c.max_d; ++d) {
    const std::vector<std::size_t> &idx = c.by_d[static_cast<std::size_t>(d)];
    if (idx.size() < 32)
      continue; // too few to time meaningfully
    // Enough repetitions that each measurement is milliseconds, not noise.
    const std::size_t work = idx.size() * std::max<std::size_t>(
                                              1, QuadPartitions(d));
    const int iters = static_cast<int>(
        std::clamp<std::size_t>(20'000'000 / std::max<std::size_t>(work, 1), 1,
                                2000));
    const auto n = static_cast<double>(idx.size()) * iters;

    const Timing hashed = TimeIt(repeats, [&] {
      for (int t = 0; t < iters; ++t)
        for (const std::size_t i : idx)
          meet_hashed<trunc_key_t>(d, c.nu.data() + c.at[i], scratch, sink);
    });
    Timing direct;
    direct.seconds = 0;
    if (d <= kMaxSweepD) {
      direct = TimeIt(repeats, [&] {
        for (int t = 0; t < iters; ++t)
          for (const std::size_t i : idx)
            DirectAt<4>(d, c.nu.data() + c.at[i], sink);
      });
    }

    const Ceiling ceil{DirectCeilingCycles(d), "hoist + load-bound table"};
    const double dcyc = direct.counters.cycles
                            ? static_cast<double>(direct.counters.cycles) / n
                            : 0.0;
    const double hcyc = hashed.counters.cycles
                            ? static_cast<double>(hashed.counters.cycles) / n
                            : 0.0;
    std::printf("  %3d %9zu", d, idx.size());
    if (direct.seconds > 0)
      std::printf(" %10.1f", direct.seconds * 1e9 / n);
    else
      std::printf(" %10s", "-");
    std::printf(" %10.1f", hashed.seconds * 1e9 / n);
    const char *win = direct.seconds <= 0 ? "hashed"
                      : direct.seconds < hashed.seconds ? "direct"
                                                        : "hashed";
    std::printf(" %8s %10.1f %10.1f", win, dcyc, hcyc);
    if (dcyc > 0)
      std::printf(" %8.1f%%", 100.0 * ceil.achieved(dcyc));
    std::printf("\n");
  }

  std::fprintf(stderr, "sink %zu\n", hits);
  return 0;
}
