// Unit tests for the per-norm E_4 enumeration.
//
// The admissible-norm counts below are the load-bearing fixtures: they were
// produced independently by the numpy sieve in drafts/verify_gpu_algos.py, so
// agreement here is a cross-language check on the generator, which is what
// every exhaustiveness claim rests on.
#include "cpu/enumerate/enumerate.h"

#include "cpu/gaussian/gaussian.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <vector>

namespace {

using pte::cpu::Gauss;
using pte::cpu::E4Sol;

struct Counts {
  std::size_t norms = 0;
  std::size_t pairs = 0;
};

Counts CountAdmissible(std::int64_t bound) {
  Counts c;
  pte::cpu::Enumerator(bound).for_each_norm([&](const pte::cpu::DecomposedNorm &n) {
    const std::size_t d = pte::cpu::point_count(n.factors);
    ++c.norms;
    c.pairs += d * (d - 1) / 2;
  });
  return c;
}

TEST(Enumerate, SplitPrimes) {
  const auto primes = pte::cpu::split_primes_up_to(100);
  const std::vector<std::uint64_t> expected = {5,  13, 17, 29, 37, 41,
                                               53, 61, 73, 89, 97};
  EXPECT_EQ(primes, expected);
}

TEST(Enumerate, PointCountMatchesTheDivisorFormula) {
  EXPECT_EQ(pte::cpu::point_count(std::vector<pte::cpu::factor_t>{5, 13, 17}), 4u);     // 8/2
  EXPECT_EQ(pte::cpu::point_count(std::vector<pte::cpu::factor_t>{5, 5, 13, 17}), 6u);  // 3*4/2
  EXPECT_EQ(pte::cpu::point_count(std::vector<pte::cpu::factor_t>{5, 5, 5, 13}), 4u);   // 4*2/2
  EXPECT_EQ(pte::cpu::point_count(std::vector<pte::cpu::factor_t>{5, 13, 17, 29}), 8u); // 16/2
  EXPECT_EQ(pte::cpu::point_count(std::vector<pte::cpu::factor_t>{5, 5}), 2u);          // ceil(3/2)
}

// Cross-language fixture: these exact numbers come from the independent numpy
// sieve in drafts/verify_gpu_algos.py.
TEST(Enumerate, AdmissibleNormCountsMatchThePythonSieve) {
  const Counts at1e6 = CountAdmissible(1000000);
  EXPECT_EQ(at1e6.norms, 11083u);
  EXPECT_EQ(at1e6.pairs, 99354u);

  const Counts at1e7 = CountAdmissible(10000000);
  EXPECT_EQ(at1e7.norms, 131437u);
  EXPECT_EQ(at1e7.pairs, 1339509u);
}

TEST(Enumerate, EveryAdmissibleNormIsEmittedOnceAndIsWellFormed) {
  std::vector<std::int64_t> seen;
  const pte::cpu::Enumerator e(200000);
  e.for_each_norm([&](const pte::cpu::DecomposedNorm &n) {
    const auto primes = e.primes_of(n);
    std::int64_t product = 1;
    for (const pte::cpu::factor_t p : primes) {
      EXPECT_EQ(p % 4, 1u) << "non-split prime in factorization";
      product *= static_cast<std::int64_t>(p);
    }
    EXPECT_EQ(product, n.l1);
    EXPECT_LE(n.l1, 200000);
    EXPECT_GE(pte::cpu::point_count(n.factors), 4u);
    EXPECT_TRUE(std::is_sorted(primes.begin(), primes.end()));
    EXPECT_TRUE(std::is_sorted(n.factors.begin(), n.factors.end()))
        << "factor indices must ascend with the primes they name";
    seen.push_back(n.l1);
  });
  const std::size_t before = seen.size();
  std::sort(seen.begin(), seen.end());
  seen.erase(std::unique(seen.begin(), seen.end()), seen.end());
  EXPECT_EQ(seen.size(), before) << "a norm was emitted more than once";
  EXPECT_EQ(seen.front(), 1105); // 5*13*17, the smallest admissible norm
}

TEST(Enumerate, NonSquarefreeNormsAreIncluded) {
  // 5^2*13*17 = 5525 must be emitted: it is exactly the case where one pair can
  // use the rational factor 5 while another uses q^2, so dropping it would
  // reopen a coverage hole.
  bool found = false;
  const pte::cpu::Enumerator e(6000);
  e.for_each_norm([&](const pte::cpu::DecomposedNorm &n) {
    if (n.l1 == 5525) {
      found = true;
      const std::vector<pte::cpu::factor_t> expected = {5, 5, 13, 17};
      EXPECT_EQ(e.primes_of(n), expected);
      EXPECT_EQ(pte::cpu::point_count(n.factors), 6u);
    }
  });
  EXPECT_TRUE(found);
}

TEST(Enumerate, MeetFindsAPlantedCollision) {
  // 1^2 + 8^2 = 4^2 + 7^2 = 65. Points 0,3 pair against 1,2.
  const std::vector<std::int64_t> nu = {1, 4, 7, 8, 11};
  const auto hits = pte::cpu::meet_oracle(nu);
  ASSERT_EQ(hits.size(), 1u);
  EXPECT_EQ(hits[0].a, 0);
  EXPECT_EQ(hits[0].b, 3);
  EXPECT_EQ(hits[0].c, 1);
  EXPECT_EQ(hits[0].d, 2);
  EXPECT_TRUE(pte::cpu::confirm_candidate(nu, hits[0]));
}

TEST(Enumerate, MeetFindsNothingWhenThereIsNothing) {
  const std::vector<std::int64_t> nu = {1, 2, 4, 8, 16, 32};
  EXPECT_TRUE(pte::cpu::meet_oracle(nu).empty());
  EXPECT_TRUE(pte::cpu::meet_oracle({}).empty());
  EXPECT_TRUE(pte::cpu::meet_oracle({3, 5, 7}).empty()); // fewer than 4
}

TEST(Enumerate, VerifyExactRejectsMalformedQuadruples) {
  const std::vector<std::int64_t> nu = {1, 4, 7, 8};
  EXPECT_TRUE(pte::cpu::confirm_candidate(nu, {0, 3, 1, 2}));
  EXPECT_FALSE(pte::cpu::confirm_candidate(nu, {0, 0, 1, 2})); // repeated index
  EXPECT_FALSE(pte::cpu::confirm_candidate(nu, {0, 3, 1, 9})); // out of range
  EXPECT_FALSE(pte::cpu::confirm_candidate(nu, {0, 1, 2, 3})); // keys differ
}

TEST(Enumerate, PrimitivityKeepsPartiallyImprimitiveQuadruples) {
  // L1 = 5^2*13*17: |U| = 6, of which exactly two points are divisible by 5.
  // Those two are 5 * U(13*17).
  const std::vector<pte::cpu::factor_t> primes = {5, 5, 13, 17};
  const auto points = pte::cpu::point_set(primes);
  ASSERT_EQ(points.size(), 6u);

  std::vector<int> imprimitive;
  for (int i = 0; i < static_cast<int>(points.size()); ++i) {
    if (points[i].re % 5 == 0 && points[i].im % 5 == 0)
      imprimitive.push_back(i);
  }
  ASSERT_EQ(imprimitive.size(), 2u) << "expected a third of the points";

  std::vector<int> primitive;
  for (int i = 0; i < static_cast<int>(points.size()); ++i) {
    if (std::find(imprimitive.begin(), imprimitive.end(), i) ==
        imprimitive.end())
      primitive.push_back(i);
  }
  ASSERT_GE(primitive.size(), 2u);

  // Two imprimitive plus two primitive points: kept (the S3 case).
  const E4Sol mixed{imprimitive[0], imprimitive[1], primitive[0],
                        primitive[1]};
  EXPECT_TRUE(pte::cpu::is_primitive(points, primes, mixed));

  // A squarefree norm can never be imprimitive.
  const std::vector<pte::cpu::factor_t> sqfree = {5, 13, 17, 29};
  const auto sqfree_points = pte::cpu::point_set(sqfree);
  EXPECT_TRUE(pte::cpu::is_primitive(sqfree_points, sqfree, {0, 1, 2, 3}));
}

TEST(Enumerate, PrimitivityRejectsFullyImprimitiveQuadruples) {
  // Build a norm whose imprimitive family is large enough to hold all four
  // points: L1 = 5^2 * 13 * 17 * 29 has |U(13*17*29)| = 4 points scaled by 5.
  const std::vector<pte::cpu::factor_t> primes = {5, 5, 13, 17, 29};
  const auto points = pte::cpu::point_set(primes);

  std::vector<int> imprimitive;
  for (int i = 0; i < static_cast<int>(points.size()); ++i) {
    if (points[i].re % 5 == 0 && points[i].im % 5 == 0)
      imprimitive.push_back(i);
  }
  ASSERT_GE(imprimitive.size(), 4u);

  const E4Sol all_imprimitive{imprimitive[0], imprimitive[1],
                                  imprimitive[2], imprimitive[3]};
  EXPECT_FALSE(pte::cpu::is_primitive(points, primes, all_imprimitive));
}

TEST(Enumerate, ExemplarNormYieldsTheKnownSolution) {
  // The E_4 exemplar of paper.md section 7, reached by targeting its norm
  // rather than sweeping to 3.3e10.
  const std::vector<pte::cpu::factor_t> primes = {13, 89, 173, 233, 709};
  const auto points = pte::cpu::point_set(primes);
  ASSERT_EQ(points.size(), 16u);

  const auto nu = pte::cpu::nu_values(points);
  const auto hits = pte::cpu::meet_oracle(nu);
  ASSERT_FALSE(hits.empty()) << "the known E_4 solution was not found";

  // The published solution has pair keys built from these nu values.
  const std::vector<std::int64_t> want = {11810488155, 32679910245, 15352988955,
                                          31172903205};
  bool matched = false;
  for (const E4Sol &q : hits) {
    std::vector<std::int64_t> got = {nu[q.a], nu[q.b], nu[q.c], nu[q.d]};
    std::vector<std::int64_t> expect = want;
    std::sort(got.begin(), got.end());
    std::sort(expect.begin(), expect.end());
    if (got == expect) {
      matched = true;
      EXPECT_TRUE(pte::cpu::confirm_candidate(nu, q));
      EXPECT_TRUE(pte::cpu::is_primitive(points, primes, q));
    }
  }
  EXPECT_TRUE(matched) << "hits found, but not the published quadruple";
}


// C11: the bound/65 factor-base cutoff is load-bearing and was untested.
//
// The argument for it is that a prime with exponent 1 needs a cofactor with
// prod(b+1) >= 4, and the smallest such cofactor is 5*13 = 65. If that were off
// by one, admissible norms would be skipped -- silently, since a search that
// misses solutions produces nothing and looks exactly like a search that found
// nothing. So this checks it by brute force rather than by re-reading the
// argument: enumerate with a *much* larger factor base and demand the same set.
TEST(Enumerate, FactorBaseCutoffLosesNothing) {
  constexpr std::int64_t kBound = 1000000;

  std::vector<std::int64_t> from_cutoff;
  pte::cpu::Enumerator(kBound).for_each_norm(
      [&](const pte::cpu::DecomposedNorm &n) { from_cutoff.push_back(n.l1); });
  std::sort(from_cutoff.begin(), from_cutoff.end());

  // Independent DFS over a factor base 13x wider than the one under test.
  const auto wide = pte::cpu::split_primes_up_to(kBound / 5);
  std::vector<std::int64_t> from_wide;
  const auto dfs = [&](auto &&self, std::size_t start, std::int64_t value,
                       std::size_t divisors) -> void {
    for (std::size_t i = start; i < wide.size(); ++i) {
      const auto p = static_cast<std::int64_t>(wide[i]);
      if (value > kBound / p)
        break;
      std::int64_t v = value * p;
      std::size_t e = 1;
      while (true) {
        if (divisors * (e + 1) >= 8)
          from_wide.push_back(v);
        self(self, i + 1, v, divisors * (e + 1));
        if (v > kBound / p)
          break;
        v *= p;
        ++e;
      }
    }
  };
  dfs(dfs, 0, 1, 1);
  std::sort(from_wide.begin(), from_wide.end());

  EXPECT_EQ(from_cutoff, from_wide)
      << "the bound/65 cutoff dropped admissible norms";
  EXPECT_EQ(from_cutoff.size(), 11083u);
}

// The cutoff is not merely correct, it is *exactly* tight: at B = 10^6 the
// largest prime factor occurring in any admissible norm is 15377, against a
// cutoff of 10^6/65 = 15384. Pinning the extremal case means a future change to
// kSmallestCofactor cannot quietly gain slack it does not have.
TEST(Enumerate, ExtremalNormUsesNearlyTheWholeFactorBase) {
  constexpr std::int64_t kBound = 1000000;
  std::uint64_t largest = 0;
  bool saw_extremal = false;
  const pte::cpu::Enumerator e(kBound);
  e.for_each_norm([&](const pte::cpu::DecomposedNorm &n) {
    largest = std::max(largest, e.prime(n.factors.back()));
    if (n.l1 == 999505)
      saw_extremal = true;
  });
  EXPECT_EQ(largest, 15377u);
  EXPECT_LE(largest, static_cast<std::uint64_t>(kBound / 65));
  EXPECT_TRUE(saw_extremal) << "5 * 13 * 15377 = 999505 must be enumerated";
}

// C6 was that Norm::primes aliased live DFS state. It still does -- the
// alternative is a heap allocation per norm in a 2e10-iteration loop -- but the
// contract is now documented and checked, so a caller who copies is safe.
TEST(Enumerate, NormIsUsableWhenCopiedInTheCallback) {
  std::vector<pte::cpu::DecomposedNorm> kept;
  const pte::cpu::Enumerator e(20000);
  e.for_each_norm([&](const pte::cpu::DecomposedNorm &n) { kept.push_back(n); });
  ASSERT_FALSE(kept.empty());
  for (const auto &n : kept) {
    std::int64_t product = 1;
    for (const pte::cpu::factor_t p : e.primes_of(n))
      product *= static_cast<std::int64_t>(p);
    EXPECT_EQ(product, n.l1) << "copied DecomposedNorm did not survive the DFS";
  }
}


TEST(Enumerate, E4SolCanonicalFormIsStableAndOrderInsensitive) {
  const E4Sol want{0, 3, 1, 2};
  // Every way of writing the same solution -- swap within a pair, swap the
  // pairs, both -- must land on one representative.
  for (const E4Sol v : {E4Sol{0, 3, 1, 2}, E4Sol{3, 0, 1, 2}, E4Sol{0, 3, 2, 1},
                        E4Sol{1, 2, 0, 3}, E4Sol{2, 1, 3, 0}})
    EXPECT_EQ(v.canonical(), want) << "not canonical: " << v.a << v.b << v.c
                                   << v.d;
  EXPECT_EQ(want.canonical().canonical(), want) << "canonical is idempotent";

  // Different pairings of the same 4-subset are different solutions and must
  // NOT be merged: a 4-subset can satisfy more than one pairing.
  EXPECT_NE(E4Sol({0, 1, 2, 3}).canonical(), E4Sol({0, 2, 1, 3}).canonical());
  EXPECT_NE(E4Sol({0, 2, 1, 3}).canonical(), E4Sol({0, 3, 1, 2}).canonical());
}

} // namespace