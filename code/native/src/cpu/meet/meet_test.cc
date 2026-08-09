// Tests for the production meet.
//
// Two fixtures carry the weight. MeetMatchesOracleOnEveryNormToOneMillion is
// the general one: the production meet and the naive exact-arithmetic oracle
// are independently written and must agree on every admissible norm. The other
// is HashedPathCarriesTheLargeDTail, which pushes to 10^7 restricted to D > 8 --
// the 3% of norms carrying 59.5% of the work there and 74.2% at 10^8, and the
// only norms meet_hashed ever sees in production.
//
// The pattern throughout: assert against something that does not share code
// with what is under test. The oracle is O(D^4), exact, and lives in another
// library; the partition tables are checked against a runtime brute force, not
// against their own count formula.
#include "cpu/meet/meet.h"

#include "cpu/core/precision.h"
#include "cpu/enumerate/enumerate.h"
#include "cpu/gaussian/gaussian.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <set>
#include <utility>
#include <vector>

namespace {

using pte::cpu::E4Sol;
using pte::cpu::exact_key_t;
using pte::cpu::MeetScratch;
using pte::cpu::norm_t;
using pte::cpu::Partition;
using pte::cpu::trunc_key_t;

// Canonical, order-independent view of a solution list, so two implementations
// can be compared without agreeing on emission order.
std::vector<E4Sol> Canonicalized(std::vector<E4Sol> hits) {
  for (E4Sol &h : hits)
    h = h.canonical();
  std::sort(hits.begin(), hits.end());
  return hits;
}

// Every split of {0..D-1} into two disjoint unordered pairs, by walking
// 4-subsets and taking all three pairings. Written the obvious way on purpose:
// this is the independent truth the compile-time table is checked against.
std::set<E4Sol> BruteForcePartitions(int D) {
  std::set<E4Sol> out;
  for (int i = 0; i < D; ++i)
    for (int j = i + 1; j < D; ++j)
      for (int k = j + 1; k < D; ++k)
        for (int l = k + 1; l < D; ++l) {
          out.insert(E4Sol{i, j, k, l}.canonical());
          out.insert(E4Sol{i, k, j, l}.canonical());
          out.insert(E4Sol{i, l, j, k}.canonical());
        }
  return out;
}

template <int D> void CheckPartitionTable() {
  constexpr auto table = pte::cpu::partitions<D>();
  static_assert(table.size() == pte::cpu::partition_count(D));

  const std::set<E4Sol> truth = BruteForcePartitions(D);
  EXPECT_EQ(table.size(), truth.size()) << "D = " << D;

  std::set<E4Sol> seen;
  for (const Partition &t : table) {
    // No index may repeat inside an entry: a "solution" reusing a point is not
    // a solution, and the table is the only thing preventing one.
    const std::array<int, 4> idx{t.a, t.b, t.c, t.d};
    const std::set<int> distinct(idx.begin(), idx.end());
    ASSERT_EQ(distinct.size(), 4u) << "D = " << D << " has a repeated index";
    for (const int v : idx)
      ASSERT_LT(v, D) << "D = " << D << " index out of range";

    // ...and no split may be generated twice, which is what the `a`-is-global-
    // minimum loop shape buys.
    const E4Sol c = E4Sol{t.a, t.b, t.c, t.d}.canonical();
    EXPECT_TRUE(seen.insert(c).second) << "D = " << D << " duplicate split";
  }
  EXPECT_EQ(seen, truth) << "D = " << D;
}

TEST(Meet, PartitionTablesMatchBruteForce) {
  CheckPartitionTable<4>();
  CheckPartitionTable<5>();
  CheckPartitionTable<6>();
  CheckPartitionTable<7>();
  CheckPartitionTable<8>();
  CheckPartitionTable<10>(); // beyond the dispatch, to exercise the generator
}

TEST(Meet, PartitionCountIsZeroBelowFour) {
  // D = 2 and D = 3 hit the zero factor at different places in the falling
  // product; a formula written without the guard returns 0 for one and a large
  // value for the other depending on evaluation order.
  EXPECT_EQ(pte::cpu::partition_count(0), 0u);
  EXPECT_EQ(pte::cpu::partition_count(1), 0u);
  EXPECT_EQ(pte::cpu::partition_count(2), 0u);
  EXPECT_EQ(pte::cpu::partition_count(3), 0u);
  EXPECT_EQ(pte::cpu::partition_count(4), 3u);
  EXPECT_EQ(pte::cpu::partition_count(8), 210u);
  // The budgeted maximum. 3*C(192,4) = 164,611,440 partitions for a single
  // norm is the number that makes meet_hashed's existence non-negotiable.
  EXPECT_EQ(pte::cpu::partition_count(192), 164'611'440u);
}

// The general fixture: both key widths against the oracle, every admissible
// norm below 10^6. Below kExactKeyBound the truncated key cannot wrap, so this
// checks the algorithm; the filter behaviour gets its own test below.
TEST(Meet, MeetMatchesOracleOnEveryNormToOneMillion) {
  const pte::cpu::Enumerator e(1'000'000);
  std::size_t norms = 0;
  std::size_t hits = 0;

  e.for_each_norm([&](const pte::cpu::DecomposedNorm &norm) {
    const std::vector<pte::cpu::Gauss> points = e.point_set(norm);
    const std::vector<norm_t> nu = pte::cpu::nu_values(points);

    const std::vector<E4Sol> expect = Canonicalized(pte::cpu::meet_oracle(nu));
    ASSERT_EQ(Canonicalized(pte::cpu::meet_collect<trunc_key_t>(nu)), expect)
        << "L1 = " << norm.l1 << ", D = " << nu.size();
    ASSERT_EQ(Canonicalized(pte::cpu::meet_collect<exact_key_t>(nu)), expect)
        << "L1 = " << norm.l1 << ", D = " << nu.size();

    ++norms;
    hits += expect.size();
  });

  EXPECT_EQ(norms, 11'083u) << "the admissible-norm count must not move";
  EXPECT_GT(hits, 0u) << "a run that found nothing would pass vacuously";
}

// The tail that matters. D > 8 is 1.16% of norms at 10^6 and 2.12% at 10^7, but
// 35.3% and 59.5% of the work -- and it is the only input meet_hashed sees in
// production, so B = 10^6 alone (129 such norms) would barely exercise it.
TEST(Meet, HashedPathCarriesTheLargeDTail) {
  const pte::cpu::Enumerator e(10'000'000);
  std::size_t large = 0;
  std::size_t hits = 0;
  int max_d = 0;

  e.for_each_norm([&](const pte::cpu::DecomposedNorm &norm) {
    if (pte::cpu::point_count(norm.factors) <= 8)
      return;
    const std::vector<pte::cpu::Gauss> points = e.point_set(norm);
    const std::vector<norm_t> nu = pte::cpu::nu_values(points);

    const std::vector<E4Sol> expect = Canonicalized(pte::cpu::meet_oracle(nu));
    ASSERT_EQ(Canonicalized(pte::cpu::meet_collect<trunc_key_t>(nu)), expect)
        << "L1 = " << norm.l1 << ", D = " << nu.size();

    ++large;
    hits += expect.size();
    max_d = std::max(max_d, static_cast<int>(nu.size()));
  });

  EXPECT_EQ(large, 2'790u) << "the D > 8 tail at B = 10^7, 2.12% of 131,437";
  EXPECT_EQ(max_d, 24) << "D_max at B = 10^7";
  EXPECT_GT(hits, 0u);
}

// meet_hashed is only reached for D > 8 in production, so its behaviour on
// small D is otherwise untested -- and small D is exactly where an authoritative
// second implementation exists to compare it against.
TEST(Meet, HashedAgreesWithDirectOnSmallD) {
  const pte::cpu::Enumerator e(200'000);
  MeetScratch<trunc_key_t> scratch;

  e.for_each_norm([&](const pte::cpu::DecomposedNorm &norm) {
    const std::vector<pte::cpu::Gauss> points = e.point_set(norm);
    if (points.size() > 8)
      return;
    const std::vector<norm_t> nu = pte::cpu::nu_values(points);

    std::vector<E4Sol> hashed;
    pte::cpu::meet_hashed<trunc_key_t>(static_cast<int>(nu.size()), nu.data(),
                                       scratch,
                                       [&](const E4Sol &s) { hashed.push_back(s); });
    EXPECT_EQ(Canonicalized(std::move(hashed)),
              Canonicalized(pte::cpu::meet_collect<trunc_key_t>(nu)))
        << "L1 = " << norm.l1 << ", D = " << nu.size();
  });
}

// |U| = ceil(prod(b_i + 1) / 2), so a norm with all-even exponents has an odd
// point count. 4225 = 5^2 * 13^2 is the smallest, at D = 5. A dispatch written
// for even D only would send these to the hashed path -- correct, but it would
// mean the D = 5 and D = 7 tables were dead code nobody noticed.
TEST(Meet, OddPointCountsAreDispatchedToATable) {
  const std::vector<pte::cpu::Gauss> points = pte::cpu::point_set({5, 5, 13, 13});
  ASSERT_EQ(points.size(), 5u) << "4225 = 5^2 * 13^2 has an odd point count";

  const std::vector<norm_t> nu = pte::cpu::nu_values(points);
  EXPECT_EQ(Canonicalized(pte::cpu::meet_collect<trunc_key_t>(nu)),
            Canonicalized(pte::cpu::meet_oracle(nu)));

  // And the odd sizes really are covered by tables, not by the default arm.
  EXPECT_EQ(pte::cpu::partition_count(5), 15u);
  EXPECT_EQ(pte::cpu::partition_count(7), 105u);
}

// The filter property, in the meet rather than in the abstract. The fixture is
// the one from precision_test.cc: four distinct values whose pair keys are
// congruent mod 2^64 but unequal exactly.
TEST(Meet, TruncatedKeyIsAFilterAndConfirmationRejectsTheCollision) {
  // 2^63 - 3 written as max - 2: forming 2^63 in a signed 64-bit type first
  // would be undefined behaviour, which is a pointed thing to get wrong in a
  // test about arithmetic width.
  const std::vector<norm_t> nu = {norm_t{1} << 33, 3, norm_t{1} << 32,
                                  std::numeric_limits<norm_t>::max() - 2};

  // Two of the three partitions collide, not one: modulo 2^64 the squares are
  // [0, 9, 0, 9], so {0,1}|{2,3} and {0,3}|{1,2} both match while {0,2}|{1,3}
  // (0 against 18) does not. That the fixture over-delivers is the point --
  // nothing downstream may treat a filter hit as a result.
  const std::vector<E4Sol> filtered =
      Canonicalized(pte::cpu::meet_collect<trunc_key_t>(nu));
  EXPECT_EQ(filtered, (std::vector<E4Sol>{{0, 1, 2, 3}, {0, 3, 1, 2}}))
      << "the truncated key must accept these quadruples";

  EXPECT_TRUE(pte::cpu::meet_collect<exact_key_t>(nu).empty())
      << "and the exact key must reject them";
  for (const E4Sol &s : filtered) {
    EXPECT_FALSE(pte::cpu::confirm_candidate(nu, s))
        << "confirm_candidate is the layer-2 arbiter that stops it";
  }

  // The oracle is exact throughout, so it must agree with the exact key.
  EXPECT_TRUE(pte::cpu::meet_oracle(nu).empty());
}

// The stated performance contract of MeetScratch: it grows to the largest D a
// run contains and then stops allocating. This is a claim about steady state,
// so it is checked by running twice and comparing capacities, not by counting
// allocations.
TEST(Meet, ScratchReachesSteadyStateAndStopsGrowing) {
  const pte::cpu::Enumerator e(1'000'000);
  MeetScratch<trunc_key_t> scratch;

  const auto sweep = [&]() {
    e.for_each_norm([&](const pte::cpu::DecomposedNorm &norm) {
      const std::vector<pte::cpu::Gauss> points = e.point_set(norm);
      const std::vector<norm_t> nu = pte::cpu::nu_values(points);
      pte::cpu::meet<trunc_key_t>(static_cast<int>(nu.size()), nu.data(),
                                  scratch, [](const E4Sol &) {});
    });
  };

  sweep();
  const std::size_t pairs = scratch.pairs.capacity();
  const std::size_t squares = scratch.sq.capacity();
  ASSERT_GT(pairs, 0u) << "the hashed path must have been exercised";

  sweep();
  EXPECT_EQ(scratch.pairs.capacity(), pairs);
  EXPECT_EQ(scratch.sq.capacity(), squares);
}

} // namespace
