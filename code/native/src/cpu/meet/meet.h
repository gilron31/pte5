// The production meet: given the nu values of one norm's point set, find every
// {a,b} | {c,d} with nu_a^2 + nu_b^2 == nu_c^2 + nu_d^2.
//
// This is the inner loop of the whole search. `meet_oracle` in cpu/enumerate
// answers the same question and stays there, deliberately naive and
// deliberately sharing no code with this file -- the cross-check between them
// is worthless if they can agree by sharing a bug.
//
// ---------------------------------------------------------------------------
// Why two algorithms
// ---------------------------------------------------------------------------
//
// The naive form is O(D^4): 3*C(D,4) pair-partitions per norm. Measured over
// the full enumeration (drafts/implementation_design.md 2), that cost
// concentrates in a tail that the norm counts completely hide:
//
//   B      norms       quad-partitions   D>8 share of norms   of the *work*
//   10^6      11,083          347,919           1.16%            35.3%
//   10^7     131,437        8,295,267           2.12%            59.5%
//   10^8   1,464,114      169,730,031           2.97%            74.2%
//
// Quad-partitions grow at 20.5x per decade against 11.1x for norms, so the tail
// is not a tail for long: by the design bound essentially all the work is in
// 3% of the norms. Two algorithms, dispatched on D:
//
//   D <= 8   meet_direct<D>  -- the partition table is a compile-time constant
//                              (3 / 15 / 45 / 105 / 210 entries), so the loop
//                              is a fixed sequence of one add and one compare.
//   D >  8   meet_hashed     -- O(D^2): build the D(D-1)/2 pair keys, sort, and
//                              read off equal-key groups. At D = 32 that is
//                              107,880 partitions replaced by 496 keys.
//
// Measured over the whole corpus (every admissible norm, nu values precomputed,
// single thread, -O2 -march=native), against two baselines that separate what
// each idea bought:
//
//   B      oracle (exact,    A8: hoisted sq,   A8+C9: this file   hoist  bucket
//          no hoist, O(D^4))  trunc key, O(D^4)
//   10^6      43.9 ns/norm      32.0 ns/norm       20.4 ns/norm    1.37x  1.57x
//   10^7      75.0             55.3               34.2             1.36x  1.62x
//   10^8     123.9            100.3               48.2             1.24x  2.08x
//
// Hoisting is a flat ~1.3x; bucketing is the one that scales, tracking the D > 8
// share of the work (35.3% -> 59.5% -> 74.2%) and near its Amdahl ceiling at
// each bound. The design doc projects ~20x at B = 10^12; that is consistent
// with this trend but is not measured here.
//
// The D <= 8 cut is itself measured, not assumed. Sweeping the threshold over
// the B = 10^7 corpus: all-hashed 51.4, direct to D<=4 36.6, to D<=6 34.6, to
// D<=7 35.1, **to D<=8 33.7**, to D<=10 34.4 ns/norm. Note that a microbenchmark
// repeating a single D says the opposite -- it has hashed winning at D = 8 by
// 27% -- because holding D fixed makes the sort's branches perfectly predictable
// and never resizes the scratch. The mixed corpus is the workload; the
// microbenchmark is not.
//
// ---------------------------------------------------------------------------
// The property both rest on
// ---------------------------------------------------------------------------
//
// **nu is injective on a canonized point set**, because N(z) = m^2 + n^2 = L1
// and nu(z) = m^2 - n^2 together determine (m, n), and canonization fixes the
// signs. Two consequences, and the second is what makes meet_hashed possible at
// all:
//
//   1. Distinct pairs with equal key are automatically *disjoint*. If {a,b} and
//      {c,d} share a point then their keys agree only if nu_b = nu_d, hence
//      b = d. So no disjointness test is needed on the candidate path.
//   2. Equal-key groups are therefore exactly the solutions, with no filtering.
//
// This is a precondition, not a theorem about arbitrary input: it holds for
// point sets that have been through Gauss::canonize and deduplicated, which is
// what point_set_from returns. See meet_hashed for what happens if a caller
// violates it.
//
// ---------------------------------------------------------------------------
// D is not always even
// ---------------------------------------------------------------------------
//
// |U| = ceil(prod(b_i + 1) / 2), so every norm with all-even exponents has an
// odd point count -- 4225 = 5^2 * 13^2 is the smallest, at D = 5. Odd D is
// 1.89% / 1.24% / 0.87% of norms at B = 10^6 / 10^7 / 10^8, and it is easy to
// lose: an earlier measurement in the design doc reported "D <= 8" as
// 97.4 / 97.0 / 96.5%, which is the share of D in {4, 6, 8}. The true D <= 8
// share is 98.84 / 97.88 / 97.03%. D = 5 and D = 7 get their own tables here
// for that reason.
#pragma once

#include "cpu/core/precision.h"
#include "cpu/core/solution.h"

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace pte::cpu {

// ---------------------------------------------------------------------------
// Compile-time partition tables
// ---------------------------------------------------------------------------

// The number of ways to split D indices into two disjoint unordered pairs:
// C(D,4) 4-subsets times 3 pairings each. Zero for D < 4, including at D = 2
// and D = 3 where the falling factorial passes through the zero factor first.
constexpr std::size_t partition_count(int D) {
  if (D < 4)
    return 0;
  return static_cast<std::size_t>(D) * (D - 1) * (D - 2) * (D - 3) / 8;
}

// One split, as {a,b} | {c,d}. uint8_t because D <= 8 wherever this is used,
// which keeps the whole D = 8 table at 840 bytes and inside L1.
struct Partition {
  std::uint8_t a, b, c, d;
};

// Every {a,b} | {c,d} over D indices, each exactly once.
//
// Uniqueness comes from the shape of the loops rather than from a filter: `a`
// runs over the *global* minimum of the four, `b` over its partner, and `c` is
// forced above `a`, so a given split is generated only with `a` as its smallest
// element and the pair containing it written first. Checked against a
// brute-force enumeration for D = 4..8 and 10 in meet_test.cc rather than
// argued for and trusted.
template <int D>
constexpr std::array<Partition, partition_count(D)> partitions() {
  std::array<Partition, partition_count(D)> out{};
  std::size_t n = 0;
  for (int a = 0; a < D; ++a)
    for (int b = a + 1; b < D; ++b)
      for (int c = a + 1; c < D; ++c) {
        if (c == b)
          continue;
        for (int d = c + 1; d < D; ++d) {
          if (d == b)
            continue;
          out[n++] = {
              static_cast<std::uint8_t>(a), static_cast<std::uint8_t>(b),
              static_cast<std::uint8_t>(c), static_cast<std::uint8_t>(d)};
        }
      }
  return out;
}

// ---------------------------------------------------------------------------
// The seams: key width and hit sink
// ---------------------------------------------------------------------------

// The two admissible key widths, and no others. trunc_key_t is the fast filter
// (one 64-bit add per partition, wrapping); exact_key_t is the arbiter. Naming
// them in a concept rather than templating on any integer is deliberate: these
// two are the pair precision.h exists to keep apart.
template <typename K>
concept MeetKey = std::same_as<K, trunc_key_t> || std::same_as<K, exact_key_t>;

// Where hits go. A plain callable rather than an object with .emit(), so a
// lambda works with no adapter; the device path wraps DeviceAppendView::push in
// one.
template <typename S>
concept HitSink = std::invocable<S &, const E4Sol &>;

// Reusable working memory for meet_hashed.
//
// Owned by the caller and reused across norms, because the campaign's inner
// loop runs this ~2e10 times and a per-norm allocation is the kind of cost that
// does not show up in a profile as itself. Growing once to the largest D a
// segment contains and then staying put is the whole intent.
template <MeetKey Key> struct MeetScratch {
  // A pair and its key. Sorting this by key is the algorithm.
  struct Entry {
    Key key;
    std::uint16_t a, b; // indices into the point set; D <= 192 by design
  };

  std::vector<Key> sq;
  std::vector<Entry> pairs;
};

// ---------------------------------------------------------------------------
// The two algorithms
// ---------------------------------------------------------------------------

// Fixed-D meet over the compile-time table. `nu` must have at least D entries.
//
// The squares are hoisted once per norm (A8 in drafts/gpu_algorithms.md), which
// is why the table loop costs one add and one compare rather than two
// multiplies and an add: at D = 4 that is 12 multiplies replaced by 4, in the
// loop that is three quarters of all work in the campaign.
template <int D, MeetKey Key, HitSink S>
void meet_direct(const norm_t *nu, S &&sink) {
  Key sq[D];
  for (int i = 0; i < D; ++i) {
    const auto v = static_cast<Key>(nu[i]);
    sq[i] = v * v;
  }

  constexpr auto table = partitions<D>();
  for (const Partition &t : table) {
    if (sq[t.a] + sq[t.b] == sq[t.c] + sq[t.d])
      sink(E4Sol{t.a, t.b, t.c, t.d});
  }
}

// O(D^2) meet: build every pair key, sort, and emit each equal-key group's
// pairs of pairs. `nu` must have at least D entries.
//
// The disjointness test the header promises is absent stays absent from the
// *candidate* path -- nothing here compares indices while scanning keys. The
// check inside the emission loop runs only on an actual key collision, which is
// rare enough to be free, and exists so that a caller who hands over a
// non-canonized point set (where nu is not injective, and equal nu values would
// pair with themselves) gets nothing rather than garbage.
template <MeetKey Key, HitSink S>
void meet_hashed(int D, const norm_t *nu, MeetScratch<Key> &scratch, S &&sink) {
  if (D < 4)
    return;

  auto &sq = scratch.sq;
  auto &pairs = scratch.pairs;
  sq.resize(static_cast<std::size_t>(D));
  for (int i = 0; i < D; ++i) {
    const auto v = static_cast<Key>(nu[i]);
    sq[static_cast<std::size_t>(i)] = v * v;
  }

  pairs.clear();
  pairs.reserve(static_cast<std::size_t>(D) * (D - 1) / 2);
  for (int i = 0; i < D; ++i)
    for (int j = i + 1; j < D; ++j) {
      pairs.push_back(
          {sq[static_cast<std::size_t>(i)] + sq[static_cast<std::size_t>(j)],
           static_cast<std::uint16_t>(i), static_cast<std::uint16_t>(j)});
    }

  std::sort(pairs.begin(), pairs.end(),
            [](const auto &x, const auto &y) { return x.key < y.key; });

  for (std::size_t i = 0; i < pairs.size();) {
    std::size_t j = i + 1;
    while (j < pairs.size() && pairs[j].key == pairs[i].key)
      ++j;
    // Groups of size 1 are the overwhelming majority and cost one comparison.
    for (std::size_t p = i; p + 1 < j; ++p)
      for (std::size_t q = p + 1; q < j; ++q) {
        const auto &u = pairs[p];
        const auto &v = pairs[q];
        if (u.a == v.a || u.a == v.b || u.b == v.a || u.b == v.b)
          continue; // only reachable if nu is not injective; see above
        sink(E4Sol{u.a, u.b, v.a, v.b});
      }
    i = j;
  }
}

// ---------------------------------------------------------------------------
// Dispatch
// ---------------------------------------------------------------------------

// The meet for one norm. `nu` must have at least D entries, canonized and
// distinct.
//
// D <= 8 covers 97.0% of norms at B = 10^8 and falls to the tables; everything
// above goes to the O(D^2) path. The cut is where the threshold sweep in the
// file header puts it.
//
// D = 5 and D = 7 have tables of their own rather than falling through to the
// default arm. D = 5 earns it (0.82% of norms at 10^7, and D = 5 is where the
// direct form is at its best); D = 7 is 0.04% and the sweep cannot resolve it
// either way. They are here because the tables cost 15 and 105 entries, and
// because a dispatch written for even D only is exactly how the design doc came
// to report the D <= 8 share as the D in {4,6,8} share.
template <MeetKey Key, HitSink S>
void meet(int D, const norm_t *nu, MeetScratch<Key> &scratch, S &&sink) {
  switch (D) {
  case 4:
    return meet_direct<4, Key>(nu, sink);
  case 5:
    return meet_direct<5, Key>(nu, sink);
  case 6:
    return meet_direct<6, Key>(nu, sink);
  case 7:
    return meet_direct<7, Key>(nu, sink);
  case 8:
    return meet_direct<8, Key>(nu, sink);
  default:
    return meet_hashed<Key>(D, nu, scratch, sink);
  }
}

// Convenience form that allocates its own scratch and collects into a vector.
//
// For tests, one-off queries and the demo app. The campaign uses the sink form
// with a scratch that outlives the norm: this one allocates twice per call,
// which is exactly what the scratch exists to avoid.
template <MeetKey Key = trunc_key_t>
std::vector<E4Sol> meet_collect(std::span<const norm_t> nu) {
  MeetScratch<Key> scratch;
  std::vector<E4Sol> hits;
  meet<Key>(static_cast<int>(nu.size()), nu.data(), scratch,
            [&hits](const E4Sol &s) { hits.push_back(s); });
  return hits;
}

} // namespace pte::cpu
