// AoS -> SoA / SIMD experiment for Spine::push. Backs
// drafts/implementation_design.md 10.1.
//
// Harvests a real corpus of (parent point set, q) pairs from the cofactor DFS
// -- the same tree in the same order as Walker::walk -- then replays the
// point-set construction in five layouts and measures each against the
// sort+unique that follows it.
//
// Not part of the build: it is a measurement, not a test, and it deliberately
// reaches past Spine's public surface to reconstruct the walk. Build it by hand
// from code/native:
//
//   g++-13 -std=c++20 -O2 -DNDEBUG -I src -o /tmp/bench_layout \
//     ../../drafts/bench_point_set_layout.cc \
//     src/cpu/generate/generate.cc src/cpu/factorbase/factorbase.cc \
//     src/cpu/gaussian/gaussian.cc
//   /tmp/bench_layout 10000000 9          # bound, repeats
//
// Repeat with -march=native for the second column of the tables in 10.1. The
// headline results: the parent point set averages 1.58 points, so there is no
// vector to fill; SoA alone is *slower* than AoS; and sort+unique is ~85% of
// push. Read 10.1 before acting on any of it.
#include "cpu/factorbase/factorbase.h"
#include "cpu/gaussian/gaussian.h"
#include "cpu/generate/generate.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <numeric>
#include <vector>

using namespace pte::cpu;
using Clock = std::chrono::steady_clock;

namespace {

// ---------------------------------------------------------------------------
// Corpus
// ---------------------------------------------------------------------------

struct Corpus {
  std::vector<Gauss> pool;      // AoS parents, packed end to end
  std::vector<std::int32_t> re; // same data, SoA
  std::vector<std::int32_t> im;
  std::vector<std::uint32_t> start;
  std::vector<std::int32_t> len;
  std::vector<Gauss> q;
  std::size_t max_len = 0;
};

// Mirrors Walker::walk, but records the inputs of every Spine::push instead of
// emitting norms. Same tree, same order, so the size distribution is the real
// one rather than a guess.
Corpus Harvest(norm_t bound, std::size_t cap) {
  Corpus c;
  Generator g(bound);
  const FactorBase &base = g.base();
  Spine spine(base);

  const auto record = [&](factor_index_t i) {
    if (c.q.size() >= cap)
      return;
    const std::span<const Gauss> p = spine.points();
    c.start.push_back(static_cast<std::uint32_t>(c.pool.size()));
    c.len.push_back(static_cast<std::int32_t>(p.size()));
    c.max_len = std::max(c.max_len, p.size());
    for (const Gauss &z : p) {
      c.pool.push_back(z);
      c.re.push_back(z.re);
      c.im.push_back(z.im);
    }
    c.q.push_back(base.point(i));
  };

  for (std::size_t segment = 0; segment < g.segment_count(); ++segment) {
    if (c.q.size() >= cap)
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
        std::size_t exponent = 1;
        while (true) {
          record(static_cast<factor_index_t>(i));
          spine.push(static_cast<factor_index_t>(i));
          self(self, i + 1, v);
          if (v > bound / p)
            break;
          v *= p;
          ++exponent;
        }
        for (std::size_t e = 0; e < exponent; ++e)
          spine.pop();
      }
    };

    norm_t v = top;
    std::size_t exponent = 1;
    while (true) {
      record(static_cast<factor_index_t>(segment));
      spine.push(static_cast<factor_index_t>(segment));
      dfs(dfs, 0, v);
      if (v > bound / top)
        break;
      v *= top;
      ++exponent;
    }
    for (std::size_t e = 0; e < exponent; ++e)
      spine.pop();
  }
  return c;
}

// ---------------------------------------------------------------------------
// Kernels: parent -> 2*parent canonized products
// ---------------------------------------------------------------------------

// V0 -- exactly what Spine::push does today: int64 intermediates, AoS,
// interleaved stores (a, b, a, b, ...).
inline void KernelV0(const Gauss *src, int n, Gauss q, Gauss *dst) {
  int m = 0;
  for (int k = 0; k < n; ++k) {
    Gauss a, b;
    src[k].mul_and_conj_mul(q, a, b);
    dst[m++] = a.canonize();
    dst[m++] = b.canonize();
  }
}

// V1 -- same arithmetic, but the two products go to two separate half-ranges.
// The result is a permutation of V0's, which a sort makes irrelevant.
inline void KernelV1(const Gauss *src, int n, Gauss q, Gauss *dst) {
  for (int k = 0; k < n; ++k) {
    Gauss a, b;
    src[k].mul_and_conj_mul(q, a, b);
    dst[k] = a.canonize();
    dst[n + k] = b.canonize();
  }
}

// V2 -- SoA, split stores, int64 intermediates.
inline void KernelV2(const std::int32_t *sre, const std::int32_t *sim, int n,
                     std::int32_t qre, std::int32_t qim, std::int32_t *dre,
                     std::int32_t *dim) {
  for (int k = 0; k < n; ++k) {
    const std::int64_t x = std::int64_t{sre[k]} * qre;
    const std::int64_t y = std::int64_t{sre[k]} * qim;
    const std::int64_t w = std::int64_t{sim[k]} * qre;
    const std::int64_t v = std::int64_t{sim[k]} * qim;
    const auto ar = static_cast<std::int32_t>(x - v);
    const auto ai = static_cast<std::int32_t>(y + w);
    const auto br = static_cast<std::int32_t>(x + v);
    const auto bi = static_cast<std::int32_t>(w - y);
    const std::int32_t aa = ar < 0 ? -ar : ar, ab = ai < 0 ? -ai : ai;
    const std::int32_t ba = br < 0 ? -br : br, bb = bi < 0 ? -bi : bi;
    dre[k] = aa >= ab ? aa : ab;
    dim[k] = aa >= ab ? ab : aa;
    dre[n + k] = ba >= bb ? ba : bb;
    dim[n + k] = ba >= bb ? bb : ba;
  }
}

// V3 -- SoA, split stores, uint32 wraparound intermediates.
//
// The correctness argument is the same one that licenses the truncated meet
// key: arithmetic mod 2^32 is a ring homomorphism, so the low 32 bits of
// (x - v) are the low 32 bits of the true product difference, and the true
// value fits in int32 (|Re(zq)| <= sqrt(N(zq)) <= 10^6). The intermediates x
// and v overflow; the answer does not, and only the answer is read.
inline void KernelV3(const std::int32_t *sre, const std::int32_t *sim, int n,
                     std::int32_t qre, std::int32_t qim, std::int32_t *dre,
                     std::int32_t *dim) {
  const auto uqre = static_cast<std::uint32_t>(qre);
  const auto uqim = static_cast<std::uint32_t>(qim);
  for (int k = 0; k < n; ++k) {
    const auto sr = static_cast<std::uint32_t>(sre[k]);
    const auto si = static_cast<std::uint32_t>(sim[k]);
    const std::uint32_t x = sr * uqre;
    const std::uint32_t y = sr * uqim;
    const std::uint32_t w = si * uqre;
    const std::uint32_t v = si * uqim;
    const auto ar = static_cast<std::int32_t>(x - v);
    const auto ai = static_cast<std::int32_t>(y + w);
    const auto br = static_cast<std::int32_t>(x + v);
    const auto bi = static_cast<std::int32_t>(w - y);
    const std::int32_t aa = ar < 0 ? -ar : ar, ab = ai < 0 ? -ai : ai;
    const std::int32_t ba = br < 0 ? -br : br, bb = bi < 0 ? -bi : bi;
    dre[k] = aa >= ab ? aa : ab;
    dim[k] = aa >= ab ? ab : aa;
    dre[n + k] = ba >= bb ? ba : bb;
    dim[n + k] = ba >= bb ? bb : ba;
  }
}

// V4 -- AoS in, packed uint64 out: (re << 32) | im, which for non-negative
// canonized coordinates orders identically to Gauss::operator<. Lets the sort
// that follows compare one word instead of two fields.
inline void KernelV4(const Gauss *src, int n, Gauss q, std::uint64_t *dst) {
  for (int k = 0; k < n; ++k) {
    Gauss a, b;
    src[k].mul_and_conj_mul(q, a, b);
    const Gauss ca = a.canonize(), cb = b.canonize();
    dst[k] = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(ca.re))
              << 32) |
             static_cast<std::uint32_t>(ca.im);
    dst[n + k] =
        (static_cast<std::uint64_t>(static_cast<std::uint32_t>(cb.re)) << 32) |
        static_cast<std::uint32_t>(cb.im);
  }
}

// ---------------------------------------------------------------------------
// Timing
// ---------------------------------------------------------------------------

struct Timing {
  double min_ns = 0;
  double med_ns = 0;
};

template <typename F> Timing Bench(int repeats, std::size_t units, F &&f) {
  std::vector<double> t;
  t.reserve(static_cast<std::size_t>(repeats));
  for (int r = 0; r < repeats; ++r) {
    const auto a = Clock::now();
    f();
    const auto b = Clock::now();
    t.push_back(std::chrono::duration<double, std::nano>(b - a).count() /
                static_cast<double>(units));
  }
  std::sort(t.begin(), t.end());
  return {t.front(), t[t.size() / 2]};
}

} // namespace

int main(int argc, char **argv) {
  const norm_t bound = argc > 1 ? std::atoll(argv[1]) : 10'000'000;
  const int repeats = argc > 2 ? std::atoi(argv[2]) : 9;

  std::fprintf(stderr, "harvesting corpus at B = %lld ...\n",
               static_cast<long long>(bound));
  const Corpus c = Harvest(bound, 4'000'000);
  const std::size_t pushes = c.q.size();
  const std::size_t points =
      std::accumulate(c.len.begin(), c.len.end(), std::size_t{0},
                      [](std::size_t s, std::int32_t n) {
                        return s + static_cast<std::size_t>(n);
                      });
  std::printf("corpus: %zu pushes, %zu parent points, mean parent D = %.2f, "
              "max %zu\n",
              pushes, points, static_cast<double>(points) / (double)pushes,
              c.max_len);

  {
    std::vector<std::size_t> hist(c.max_len + 1, 0);
    for (const std::int32_t n : c.len)
      ++hist[static_cast<std::size_t>(n)];
    std::printf("parent-size histogram (the SIMD trip count):\n");
    for (std::size_t n = 0; n < hist.size(); ++n)
      if (hist[n])
        std::printf("  D_parent = %2zu : %8zu pushes (%5.1f%%)\n", n, hist[n],
                    100.0 * (double)hist[n] / (double)pushes);
  }

  std::vector<Gauss> out(2 * c.max_len + 8);
  std::vector<std::int32_t> ore(2 * c.max_len + 8), oim(2 * c.max_len + 8);
  std::vector<std::uint64_t> opk(2 * c.max_len + 8);

  // Correctness: V1..V4 must be permutations of V0's output.
  {
    bool ok = true;
    std::vector<Gauss> a, b;
    for (std::size_t i = 0; i < pushes && ok; ++i) {
      const int n = c.len[i];
      const Gauss *src = c.pool.data() + c.start[i];
      KernelV0(src, n, c.q[i], out.data());
      a.assign(out.begin(), out.begin() + 2 * n);
      KernelV1(src, n, c.q[i], out.data());
      b.assign(out.begin(), out.begin() + 2 * n);
      KernelV3(c.re.data() + c.start[i], c.im.data() + c.start[i], n,
               c.q[i].re, c.q[i].im, ore.data(), oim.data());
      std::vector<Gauss> d;
      for (int k = 0; k < 2 * n; ++k)
        d.push_back(Gauss{ore[k], oim[k]});
      std::sort(a.begin(), a.end());
      std::sort(b.begin(), b.end());
      std::sort(d.begin(), d.end());
      ok = (a == b) && (a == d);
    }
    std::printf("V1/V3 agree with V0 on every push: %s\n", ok ? "yes" : "NO");
    if (!ok)
      return 1;
  }

  std::uint64_t sink = 0;
  const auto sweep = [&](const char *name, auto &&body) {
    const Timing t = Bench(repeats, pushes, body);
    std::printf("  %-34s %7.2f ns/push (min)  %7.2f (median)\n", name, t.min_ns,
                t.med_ns);
  };

  std::printf("\nmultiply + canonize only:\n");
  sweep("V0 AoS, int64, interleaved", [&] {
    for (std::size_t i = 0; i < pushes; ++i) {
      KernelV0(c.pool.data() + c.start[i], c.len[i], c.q[i], out.data());
      sink += static_cast<std::uint64_t>(out[0].re);
    }
  });
  sweep("V1 AoS, int64, split stores", [&] {
    for (std::size_t i = 0; i < pushes; ++i) {
      KernelV1(c.pool.data() + c.start[i], c.len[i], c.q[i], out.data());
      sink += static_cast<std::uint64_t>(out[0].re);
    }
  });
  sweep("V2 SoA, int64, split stores", [&] {
    for (std::size_t i = 0; i < pushes; ++i) {
      KernelV2(c.re.data() + c.start[i], c.im.data() + c.start[i], c.len[i],
               c.q[i].re, c.q[i].im, ore.data(), oim.data());
      sink += static_cast<std::uint64_t>(ore[0]);
    }
  });
  sweep("V3 SoA, uint32 wrap, split stores", [&] {
    for (std::size_t i = 0; i < pushes; ++i) {
      KernelV3(c.re.data() + c.start[i], c.im.data() + c.start[i], c.len[i],
               c.q[i].re, c.q[i].im, ore.data(), oim.data());
      sink += static_cast<std::uint64_t>(ore[0]);
    }
  });
  sweep("V4 AoS in, packed uint64 out", [&] {
    for (std::size_t i = 0; i < pushes; ++i) {
      KernelV4(c.pool.data() + c.start[i], c.len[i], c.q[i], opk.data());
      sink += opk[0];
    }
  });

  std::printf("\nmultiply + canonize + sort + unique (the whole of push):\n");
  sweep("V0 + std::sort<Gauss> + unique", [&] {
    for (std::size_t i = 0; i < pushes; ++i) {
      const int n = c.len[i];
      KernelV0(c.pool.data() + c.start[i], n, c.q[i], out.data());
      std::sort(out.data(), out.data() + 2 * n);
      const auto e = std::unique(out.data(), out.data() + 2 * n);
      sink += static_cast<std::uint64_t>(e - out.data());
    }
  });
  sweep("V4 + std::sort<uint64> + unique", [&] {
    for (std::size_t i = 0; i < pushes; ++i) {
      const int n = c.len[i];
      KernelV4(c.pool.data() + c.start[i], n, c.q[i], opk.data());
      std::sort(opk.data(), opk.data() + 2 * n);
      const auto e = std::unique(opk.data(), opk.data() + 2 * n);
      sink += static_cast<std::uint64_t>(e - opk.data());
    }
  });

  std::printf("\nsort + unique alone (on V0 output, already in cache):\n");
  sweep("std::sort<Gauss> + unique", [&] {
    for (std::size_t i = 0; i < pushes; ++i) {
      const int n = c.len[i];
      std::memcpy(out.data(), c.pool.data() + c.start[i],
                  sizeof(Gauss) * static_cast<std::size_t>(n));
      std::memcpy(out.data() + n, c.pool.data() + c.start[i],
                  sizeof(Gauss) * static_cast<std::size_t>(n));
      std::sort(out.data(), out.data() + 2 * n);
      const auto e = std::unique(out.data(), out.data() + 2 * n);
      sink += static_cast<std::uint64_t>(e - out.data());
    }
  });

  std::fprintf(stderr, "sink %llu\n", (unsigned long long)sink);
  return 0;
}
