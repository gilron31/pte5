// B4 -- the sort inside Spine::push, against a sorting-network ceiling.
//
// bench_walk_phases puts Spine::push at ~40 ns/norm, flat across three decades.
// drafts/implementation_design.md 10.1 then found that ~85% of that is
// `std::sort` + `std::unique` on a **three-element** array: the parent point set
// averages 1.58 points, so push sorts 2 x 1.58 products.
//
// std::sort is the wrong tool at that size and the counters say why -- it is a
// comparison sort, so it branches on data, and a branch on data at n = 3 is a
// coin flip the predictor cannot learn. A sorting network has no data-dependent
// branches at all: it is a fixed sequence of min/max pairs. This bench measures
// how much of the gap that closes, and prints the network's own floor so the
// answer is "x% of achievable" rather than "some percent faster".
//
// Deliberately does NOT change Spine. Reordering points changes the catalog's
// nu column order, so adopting any of this is a separate, gated change.
//
//   ./bench_sort_ceiling <bound> <repeats>
#include "bench/ceiling.h"
#include "bench/harness.h"
#include "cpu/core/precision.h"
#include "cpu/factorbase/factorbase.h"
#include "cpu/gaussian/gaussian.h"
#include "cpu/generate/generate.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace pte::cpu;
using pte::bench::Ceiling;
using pte::bench::Machine;
using pte::bench::Suite;

namespace {

// ---------------------------------------------------------------------------
// Corpus: the arrays push actually sorts
// ---------------------------------------------------------------------------
//
// Harvested by replaying the walk and capturing, at each push, the 2*parent
// canonized products *before* the sort. Packed as uint64 with re in the high
// word: canonized coordinates are non-negative, so that ordering is identical
// to Gauss::operator< and every variant below can share one representation.

struct Corpus {
  std::vector<std::uint64_t> vals;
  std::vector<std::size_t> at;
  std::vector<int> n;
  std::size_t total = 0;
  int max_n = 0;
};

std::uint64_t Pack(const Gauss &g) {
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(g.re)) << 32) |
         static_cast<std::uint32_t>(g.im);
}

Corpus Harvest(norm_t bound, std::size_t cap) {
  Corpus c;
  Generator g(bound);
  const FactorBase &base = g.base();
  Spine spine(base);
  std::vector<std::uint64_t> scratch;

  const auto record = [&](factor_index_t i) {
    if (c.at.size() >= cap)
      return;
    const std::span<const Gauss> p = spine.points();
    const Gauss q = base.point(i);
    scratch.clear();
    for (const Gauss &z : p) {
      Gauss a, b;
      z.mul_and_conj_mul(q, a, b);
      scratch.push_back(Pack(a.canonize()));
      scratch.push_back(Pack(b.canonize()));
    }
    c.at.push_back(c.vals.size());
    c.n.push_back(static_cast<int>(scratch.size()));
    c.max_n = std::max(c.max_n, static_cast<int>(scratch.size()));
    c.total += scratch.size();
    c.vals.insert(c.vals.end(), scratch.begin(), scratch.end());
  };

  for (std::size_t segment = 0; segment < g.segment_count(); ++segment) {
    if (c.at.size() >= cap)
      break;
    const auto top =
        static_cast<norm_t>(base.prime(static_cast<factor_index_t>(segment)));
    if (top > bound)
      break;
    const auto dfs = [&](auto &&self, std::size_t begin, norm_t value) -> void {
      for (std::size_t i = begin; i < segment; ++i) {
        const auto p =
            static_cast<norm_t>(base.prime(static_cast<factor_index_t>(i)));
        if (value > bound / p)
          break;
        norm_t v = value * p;
        std::size_t e = 1;
        while (true) {
          record(static_cast<factor_index_t>(i));
          spine.push(static_cast<factor_index_t>(i));
          self(self, i + 1, v);
          if (v > bound / p)
            break;
          v *= p;
          ++e;
        }
        for (std::size_t k = 0; k < e; ++k)
          spine.pop();
      }
    };
    norm_t v = top;
    std::size_t e = 1;
    while (true) {
      record(static_cast<factor_index_t>(segment));
      spine.push(static_cast<factor_index_t>(segment));
      dfs(dfs, 0, v);
      if (v > bound / top)
        break;
      v *= top;
      ++e;
    }
    for (std::size_t k = 0; k < e; ++k)
      spine.pop();
  }
  return c;
}

// ---------------------------------------------------------------------------
// The variants
// ---------------------------------------------------------------------------

// Branchless compare-exchange. cmov, not a jump -- this is the entire point.
inline void CE(std::uint64_t &a, std::uint64_t &b) {
  const std::uint64_t lo = a < b ? a : b;
  const std::uint64_t hi = a < b ? b : a;
  a = lo;
  b = hi;
}

// Optimal-size Batcher/Knuth networks for the sizes push actually sees.
// 90% of pushes produce n = 2 or 4 (parent D of 1 or 2).
inline void Net2(std::uint64_t *v) { CE(v[0], v[1]); }
inline void Net4(std::uint64_t *v) {
  CE(v[0], v[1]); CE(v[2], v[3]);
  CE(v[0], v[2]); CE(v[1], v[3]);
  CE(v[1], v[2]);
}
inline void Net6(std::uint64_t *v) {
  CE(v[0], v[1]); CE(v[2], v[3]); CE(v[4], v[5]);
  CE(v[0], v[2]); CE(v[3], v[5]); CE(v[1], v[4]);
  CE(v[0], v[1]); CE(v[2], v[3]); CE(v[4], v[5]);
  CE(v[1], v[2]); CE(v[3], v[4]);
  CE(v[2], v[3]);
}
inline void Net8(std::uint64_t *v) {
  CE(v[0], v[1]); CE(v[2], v[3]); CE(v[4], v[5]); CE(v[6], v[7]);
  CE(v[0], v[2]); CE(v[1], v[3]); CE(v[4], v[6]); CE(v[5], v[7]);
  CE(v[1], v[2]); CE(v[5], v[6]); CE(v[0], v[4]); CE(v[3], v[7]);
  CE(v[1], v[5]); CE(v[2], v[6]);
  CE(v[1], v[4]); CE(v[3], v[6]);
  CE(v[2], v[4]); CE(v[3], v[5]);
  CE(v[3], v[4]);
}

// Compare-exchange counts, for the ceiling.
inline int NetCE(int n) {
  switch (n) {
  case 2: return 1;
  case 4: return 5;
  case 6: return 12;
  case 8: return 19;
  default: return 0;
  }
}

inline bool NetworkSort(std::uint64_t *v, int n) {
  switch (n) {
  case 0:
  case 1: return true;
  case 2: Net2(v); return true;
  case 4: Net4(v); return true;
  case 6: Net6(v); return true;
  case 8: Net8(v); return true;
  default: return false; // sizes push produces are even; others fall back
  }
}

// Branchless insertion sort: no network, but also no data-dependent branch on
// the comparison itself. The interesting middle case.
inline void InsertionCE(std::uint64_t *v, int n) {
  for (int i = 1; i < n; ++i)
    for (int j = i; j > 0; --j)
      CE(v[j - 1], v[j]);
}

std::size_t UniqueRun(std::uint64_t *v, int n) {
  return static_cast<std::size_t>(std::unique(v, v + n) - v);
}

} // namespace

int main(int argc, char **argv) {
  const norm_t bound = argc > 1 ? std::atoll(argv[1]) : 10'000'000;
  const int repeats = argc > 2 ? std::atoi(argv[2]) : 7;

  std::fprintf(stderr, "harvesting push corpus at B = %lld ...\n",
               static_cast<long long>(bound));
  const Corpus c = Harvest(bound, 8'000'000);
  if (c.at.empty()) {
    std::printf("no pushes below %lld\n", static_cast<long long>(bound));
    return 0;
  }
  const auto pushes = static_cast<double>(c.at.size());

  // Distribution, because it is the reason the answer is what it is.
  {
    std::vector<std::size_t> hist(static_cast<std::size_t>(c.max_n) + 1, 0);
    for (const int n : c.n)
      ++hist[static_cast<std::size_t>(n)];
    std::printf("=== B4 sort corpus: what Spine::push actually sorts ===\n");
    std::size_t covered = 0;
    for (std::size_t n = 0; n < hist.size(); ++n)
      if (hist[n]) {
        std::printf("  n = %2zu : %9zu pushes (%5.1f%%)%s\n", n, hist[n],
                    100.0 * static_cast<double>(hist[n]) / pushes,
                    NetCE(static_cast<int>(n)) ? "" : "   <- no network, falls back");
        if (NetCE(static_cast<int>(n)))
          covered += hist[n];
      }
    std::printf("  networks cover %.1f%% of pushes\n",
                100.0 * static_cast<double>(covered) / pushes);
  }

  // Correctness gate: every variant must produce the same deduplicated array.
  {
    std::vector<std::uint64_t> a, b;
    for (std::size_t i = 0; i < c.at.size(); ++i) {
      const int n = c.n[i];
      a.assign(c.vals.begin() + static_cast<std::ptrdiff_t>(c.at[i]),
               c.vals.begin() + static_cast<std::ptrdiff_t>(c.at[i]) + n);
      b = a;
      std::sort(a.begin(), a.end());
      a.resize(UniqueRun(a.data(), n));
      if (!NetworkSort(b.data(), n))
        std::sort(b.begin(), b.end());
      b.resize(UniqueRun(b.data(), n));
      if (a != b) {
        std::fprintf(stderr, "NETWORK WRONG at push %zu (n = %d)\n", i, n);
        return 1;
      }
    }
    std::fprintf(stderr, "network sort agrees with std::sort on all %zu "
                         "pushes\n",
                 c.at.size());
  }

  // Ceiling: a compare-exchange is cmp + 2 cmov = 3 ALU ops, plus the load and
  // store of the array. Networks have no branches, so nothing else applies.
  double ce_total = 0;
  for (const int n : c.n)
    ce_total += NetCE(n) ? NetCE(n) : n * n / 4.0; // insertion fallback estimate
  const Ceiling ceiling{
      3.0 * ce_total / pushes / Machine::kIntAluPipes,
      "3 ALU ops per compare-exchange (cmp + 2 cmov) across 4 pipes, "
      "branch-free"};

  std::vector<std::uint64_t> buf(static_cast<std::size_t>(c.max_n) + 8);
  std::uint64_t sink = 0;
  const auto load = [&](std::size_t i) {
    const int n = c.n[i];
    std::copy(c.vals.begin() + static_cast<std::ptrdiff_t>(c.at[i]),
              c.vals.begin() + static_cast<std::ptrdiff_t>(c.at[i]) + n,
              buf.begin());
    return n;
  };

  Suite s("B4 sort + unique over the real push corpus", pushes, "push");
  s.Note("  " + std::to_string(c.total) + " values, mean n = " +
         std::to_string(static_cast<double>(c.total) / pushes).substr(0, 4));
  s.SetCeiling(ceiling);
  s.Add("std::sort + unique (shipped)", [&] {
    for (std::size_t i = 0; i < c.at.size(); ++i) {
      const int n = load(i);
      std::sort(buf.data(), buf.data() + n);
      sink += UniqueRun(buf.data(), n);
    }
  });
  s.Add("branchless insertion + unique", [&] {
    for (std::size_t i = 0; i < c.at.size(); ++i) {
      const int n = load(i);
      InsertionCE(buf.data(), n);
      sink += UniqueRun(buf.data(), n);
    }
  });
  s.Add("sorting network + unique", [&] {
    for (std::size_t i = 0; i < c.at.size(); ++i) {
      const int n = load(i);
      if (!NetworkSort(buf.data(), n))
        std::sort(buf.data(), buf.data() + n);
      sink += UniqueRun(buf.data(), n);
    }
  });
  s.Add("copy only (the floor of this loop)", [&] {
    for (std::size_t i = 0; i < c.at.size(); ++i)
      sink += static_cast<std::uint64_t>(load(i)) + buf[0];
  });
  const auto &r = s.Run(repeats);

  // The "copy only" row is this bench's own overhead, not push's: the real
  // Spine::push writes its products straight into the spine buffer and sorts in
  // place, with no copy at all. Subtracting it is what makes the comparison a
  // statement about sorting rather than about memcpy.
  const double base_s = r.back().min_s;
  const double base_c = static_cast<double>(r.back().counters.cycles);
  std::printf("\n  net of the copy floor (this is the sort itself):\n");
  std::printf("  %-34s %9s %9s %9s\n", "variant", "ns/push", "cyc/push",
              "%ceil");
  for (std::size_t i = 0; i + 1 < r.size(); ++i) {
    const double dt = r[i].min_s - base_s;
    const double dc = (static_cast<double>(r[i].counters.cycles) - base_c) /
                      pushes;
    std::printf("  %-34s %9.2f %9.2f %8.1f%%\n", r[i].name.c_str(),
                dt * 1e9 / pushes, dc, 100.0 * ceiling.achieved(dc));
  }
  std::printf("\n  Note the branch-miss column above: std::sort takes ~0.2 "
              "misses/push here,\n"
              "  so at n = 3 it is NOT mispredict-bound and a branch-free "
              "network cannot win\n"
              "  by the margin that assumption would predict. Measure before "
              "rewriting Spine.\n");

  std::fprintf(stderr, "sink %llu\n", (unsigned long long)(sink + s.sink()));
  return 0;
}
