#include "cpu/tower/tower.h"
#include "cpu/tower/tower_oracle.h"

#include <gtest/gtest.h>

#include <cstdint>
#include <vector>

namespace pte::cpu {
namespace {

using U64s = std::vector<std::uint64_t>;

TowerCounts RunDp(std::uint32_t p, std::uint32_t k, std::uint32_t depth,
                std::uint32_t cap) {
  TowerParams params;
  params.p = p;
  params.k = k;
  params.depth = depth;
  params.max_valuation = cap;
  return tower_counts(params);
}

// ---------------------------------------------------------------------------
// The two rejected models are the two extreme caps at k = 1
// ---------------------------------------------------------------------------
//
// Every number below is what drafts/verify_gpu_plan.py prints for the model it
// names, so this is the assertion that the valuation cap *generalises* both
// rather than replacing them: cap 0 is "t != 0 mod p" (distinct roots) and cap
// k = 1 is no restriction at all (with multiplicity). If this test fails, the
// cap is not the right knob.

struct KOneRow {
  std::uint32_t p;
  U64s towers_cap0, states_cap0, towers_uncapped, states_uncapped;
};

TEST(TowerKOne, ReproducesBothRejectedModels) {
  const std::vector<KOneRow> rows = {
      {3, {1, 0, 0, 0, 0}, {1, 3, 1, 1, 1}, {2, 3, 3, 3, 3}, {1, 3, 3, 4, 4}},
      {5, {2, 1, 0, 0, 0}, {1, 5, 6, 1, 1}, {3, 6, 7, 7, 7}, {1, 5, 10, 6, 6}},
      {7, {3, 3, 0, 0, 0}, {1, 7, 8, 1, 1}, {4, 10, 13, 13, 13},
       {1, 7, 28, 8, 8}},
      {11, {5, 10, 0, 0, 0}, {1, 11, 12, 1, 1}, {6, 21, 31, 31, 31},
       {1, 11, 66, 12, 12}},
      {13, {6, 15, 3, 0, 0}, {1, 13, 26, 14, 1}, {7, 28, 46, 46, 46},
       {1, 13, 52, 14, 14}},
  };
  for (const KOneRow &row : rows) {
    const TowerCounts distinct = RunDp(row.p, 1, 5, 0);
    EXPECT_EQ(distinct.towers, row.towers_cap0) << "p = " << row.p;
    EXPECT_EQ(distinct.states, row.states_cap0) << "p = " << row.p;

    const TowerCounts multiplicity = RunDp(row.p, 1, 5, kNoValuationCap);
    EXPECT_EQ(multiplicity.towers, row.towers_uncapped) << "p = " << row.p;
    EXPECT_EQ(multiplicity.states, row.states_uncapped) << "p = " << row.p;
  }
}

TEST(TowerKOne, ReproducesTheLandmineCounts) {
  // The four primes <= 113 at which the distinct-root model does *not* give
  // zero depth-5 towers -- the numbers that make "zero at 25 of 29 primes"
  // a statement about the model rather than about E_5.
  const std::vector<std::pair<std::uint32_t, std::uint64_t>> expected = {
      {79, 78}, {89, 11}, {97, 207}, {113, 371}};
  for (const auto &[p, towers] : expected)
    EXPECT_EQ(RunDp(p, 1, 5, 0).towers.back(), towers) << "p = " << p;
}

// ---------------------------------------------------------------------------
// The set DP against the tree definition
// ---------------------------------------------------------------------------

TEST(TowerOracleAgreement, MatchesBruteForceAtEveryCap) {
  struct Case {
    std::uint32_t p, k, depth;
  };
  // Kept to m^depth <= ~2e4 tuples: this enumerates every tuple and walks the
  // whole tree for each, so it buys exactness with time.
  const std::vector<Case> cases = {
      {2, 1, 4}, {2, 2, 3}, {2, 3, 3}, {3, 1, 4}, {3, 2, 3},
      {3, 3, 3}, {5, 1, 4}, {5, 2, 3}, {7, 1, 3}, {7, 2, 2},
  };
  for (const Case &c : cases)
    for (std::uint32_t cap = 0; cap <= c.k; ++cap) {
      const std::uint64_t dp = RunDp(c.p, c.k, c.depth, cap).towers.back();
      const std::uint64_t brute = tower_count_brute(c.p, c.k, c.depth, cap);
      EXPECT_EQ(dp, brute) << "p = " << c.p << ", k = " << c.k
                           << ", depth = " << c.depth << ", cap = " << cap;
    }
}

// ---------------------------------------------------------------------------
// T4: the exemplar, which is what makes a null result mean anything
// ---------------------------------------------------------------------------

TEST(TowerExemplar, RecoversItsNodeValuations) {
  // ht_p of the section-7 exemplar, computed independently over Z in
  // drafts/probe_local_model.py from its 2^4 roots. The oracle sees only the
  // reduced tuple, so agreeing here means the reduction really does carry the
  // valuation -- and it is why "no depth-4 tower mod 3^5" is arithmetic rather
  // than a bug: 5 < 6.
  struct Case {
    std::uint32_t p, k, expected_floor;
  };
  const std::vector<Case> cases = {
      {7, 1, 1}, {7, 2, 2}, {7, 3, 2}, {7, 4, 2},  // ht_7 = 2
      {3, 4, 4}, {3, 5, 5}, {3, 6, 6}, {3, 7, 6},  // ht_3 = 6
      {5, 3, 3}, {5, 4, 4},                        // ht_5 = 6, out of reach here
  };
  for (const Case &c : cases) {
    const TowerOracle oracle(c.p, c.k);
    const std::vector<std::uint32_t> ells = exemplar_e4_mod(oracle.modulus());
    EXPECT_EQ(oracle.tuple_valuation_floor(ells), c.expected_floor)
        << "p = " << c.p << ", k = " << c.k;
  }
}

TEST(TowerExemplar, BoundsTheValuationFloor) {
  // A real solution reduces to a tower, so V_p(k,4) can never exceed the
  // exemplar's ht_p. At p = 7, k = 3 that is a real constraint (2 < 3) rather
  // than the vacuous V <= k.
  const ValuationFloor v = valuation_floor(7, 3, 4);
  EXPECT_LE(v.floor, 2u);
  EXPECT_TRUE(v.nondegenerate());
}

// ---------------------------------------------------------------------------
// Cheap invariants that catch most cap-handling bugs
// ---------------------------------------------------------------------------

TEST(TowerInvariants, TowersAreMonotoneInTheCap) {
  for (const std::uint32_t p : {2u, 3u, 5u})
    for (const std::uint32_t k : {2u, 3u}) {
      std::uint64_t previous = 0;
      for (std::uint32_t cap = 0; cap <= k; ++cap) {
        const std::uint64_t towers = RunDp(p, k, 4, cap).towers.back();
        EXPECT_GE(towers, previous) << "p = " << p << ", k = " << k
                                    << ", cap = " << cap;
        previous = towers;
      }
    }
}

TEST(TowerInvariants, ValuationFloorIsMonotoneInKAndDepth) {
  for (const std::uint32_t p : {3u, 5u}) {
    std::uint32_t previous_depth_floor = 0;
    for (const std::uint32_t depth : {2u, 3u, 4u}) {
      const std::uint32_t f = valuation_floor(p, 3, depth).floor;
      EXPECT_GE(f, previous_depth_floor) << "p = " << p << ", depth = " << depth;
      previous_depth_floor = f;
    }
    std::uint32_t previous_k_floor = 0;
    for (const std::uint32_t k : {1u, 2u, 3u}) {
      const std::uint32_t f = valuation_floor(p, k, 3).floor;
      EXPECT_GE(f, previous_k_floor) << "p = " << p << ", k = " << k;
      previous_k_floor = f;
    }
  }
}

// The structural claim the header rests on, and the one that decides whether a
// modulus like 5^7 is reachable at all: every state at level >= 1 has the form
// S - L, so each level's state set is a union of complete translation orbits.
// `orbit_misses` counts rotations that are not themselves states; a nonzero
// entry falsifies that, and the orbit count is what the state count would
// collapse to if a future run stored representatives instead of states.
TEST(TowerInvariants, StatesAreClosedUnderTranslation) {
  for (const std::uint32_t p : {3u, 5u})
    for (const std::uint32_t k : {1u, 2u}) {
      TowerParams params;
      params.p = p;
      params.k = k;
      params.depth = 4;
      params.count_orbits = true;
      const TowerCounts counts = tower_counts(params);
      for (std::size_t level = 0; level < counts.orbits.size(); ++level) {
        EXPECT_EQ(counts.orbit_misses[level], 0u)
            << "p = " << p << ", k = " << k << ", level = " << level;
        EXPECT_LE(counts.orbits[level], counts.states[level]);
      }
    }
}

TEST(TowerInvariants, RejectsBadParameters) {
  TowerParams params;
  params.p = 4; // not prime
  EXPECT_THROW(tower_counts(params), std::invalid_argument);
  params.p = 3;
  params.k = 0;
  EXPECT_THROW(tower_counts(params), std::invalid_argument);
  params.k = 1;
  params.depth = 0;
  EXPECT_THROW(tower_counts(params), std::invalid_argument);

  // Exact counting past 2^64 is refused, not wrapped. 3^9 = 19683 and
  // 19683^5 = 3.0e21; the guard is on m^depth, so it fires before any work.
  params.depth = 5;
  params.k = 9;
  params.p = 3;
  EXPECT_THROW(tower_counts(params), TowerTooLarge);

  // The obstruction branch does not count, so the same modulus is legal there
  // -- what bounds it is the state table instead.
  params.k = 5;
  params.depth = 3;
  params.count_towers = false;
  params.max_states = 100;
  EXPECT_THROW(tower_counts(params), TowerTooLarge);
}

TEST(TowerOracleUnits, ReducesTheExemplarDecimals) {
  EXPECT_EQ(reduce_decimal("33065996717", 1000u), 717u);
  EXPECT_EQ(reduce_decimal(kExemplarE4[0], 7u), 33065996717ull % 7u);
  EXPECT_THROW(reduce_decimal("12x4", 10u), std::invalid_argument);
  // L_4 is a perfect square, so it is a square mod every m.
  const std::uint32_t m = 2187;
  const std::uint32_t l4 = exemplar_e4_mod(m).back();
  bool is_square = false;
  for (std::uint32_t t = 0; t < m; ++t)
    if ((static_cast<std::uint64_t>(t) * t) % m == l4)
      is_square = true;
  EXPECT_TRUE(is_square);
}

} // namespace
} // namespace pte::cpu
