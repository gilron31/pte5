// Unit tests for the per-norm E_4 enumeration.
//
// The admissible-norm counts below are the load-bearing fixtures: they were
// produced independently by the numpy sieve in drafts/verify_gpu_algos.py, so
// agreement here is a cross-language check on the generator, which is what
// every exhaustiveness claim rests on.
#include "enumerate.h"

#include "gaussian.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <vector>

namespace {

using pte::cpu::Gauss;
using pte::cpu::Quadruple;

struct Counts {
  std::size_t norms = 0;
  std::size_t pairs = 0;
};

Counts CountAdmissible(std::int64_t bound) {
  Counts c;
  pte::cpu::for_each_admissible_norm(bound, [&](const pte::cpu::Norm &n) {
    const std::size_t d = pte::cpu::point_count(n.primes);
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
  EXPECT_EQ(pte::cpu::point_count({5, 13, 17}), 4u);     // 8/2
  EXPECT_EQ(pte::cpu::point_count({5, 5, 13, 17}), 6u);  // 3*4/2
  EXPECT_EQ(pte::cpu::point_count({5, 5, 5, 13}), 4u);   // 4*2/2
  EXPECT_EQ(pte::cpu::point_count({5, 13, 17, 29}), 8u); // 16/2
  EXPECT_EQ(pte::cpu::point_count({5, 5}), 2u);          // ceil(3/2)
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
  pte::cpu::for_each_admissible_norm(200000, [&](const pte::cpu::Norm &n) {
    std::int64_t product = 1;
    for (const std::uint64_t p : n.primes) {
      EXPECT_EQ(p % 4, 1u) << "non-split prime in factorization";
      product *= static_cast<std::int64_t>(p);
    }
    EXPECT_EQ(product, n.l1);
    EXPECT_LE(n.l1, 200000);
    EXPECT_GE(pte::cpu::point_count(n.primes), 4u);
    EXPECT_TRUE(std::is_sorted(n.primes.begin(), n.primes.end()));
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
  pte::cpu::for_each_admissible_norm(6000, [&](const pte::cpu::Norm &n) {
    if (n.l1 == 5525) {
      found = true;
      const std::vector<std::uint64_t> expected = {5, 5, 13, 17};
      EXPECT_EQ(n.primes, expected);
      EXPECT_EQ(pte::cpu::point_count(n.primes), 6u);
    }
  });
  EXPECT_TRUE(found);
}

TEST(Enumerate, MeetFindsAPlantedCollision) {
  // 1^2 + 8^2 = 4^2 + 7^2 = 65. Points 0,3 pair against 1,2.
  const std::vector<std::int64_t> nu = {1, 4, 7, 8, 11};
  const auto hits = pte::cpu::meet_reference(nu);
  ASSERT_EQ(hits.size(), 1u);
  EXPECT_EQ(hits[0].a, 0);
  EXPECT_EQ(hits[0].b, 3);
  EXPECT_EQ(hits[0].c, 1);
  EXPECT_EQ(hits[0].d, 2);
  EXPECT_TRUE(pte::cpu::verify_exact(nu, hits[0]));
}

TEST(Enumerate, MeetFindsNothingWhenThereIsNothing) {
  const std::vector<std::int64_t> nu = {1, 2, 4, 8, 16, 32};
  EXPECT_TRUE(pte::cpu::meet_reference(nu).empty());
  EXPECT_TRUE(pte::cpu::meet_reference({}).empty());
  EXPECT_TRUE(pte::cpu::meet_reference({3, 5, 7}).empty()); // fewer than 4
}

TEST(Enumerate, VerifyExactRejectsMalformedQuadruples) {
  const std::vector<std::int64_t> nu = {1, 4, 7, 8};
  EXPECT_TRUE(pte::cpu::verify_exact(nu, {0, 3, 1, 2}));
  EXPECT_FALSE(pte::cpu::verify_exact(nu, {0, 0, 1, 2})); // repeated index
  EXPECT_FALSE(pte::cpu::verify_exact(nu, {0, 3, 1, 9})); // out of range
  EXPECT_FALSE(pte::cpu::verify_exact(nu, {0, 1, 2, 3})); // keys differ
}

TEST(Enumerate, PrimitivityKeepsPartiallyImprimitiveQuadruples) {
  // L1 = 5^2*13*17: |U| = 6, of which exactly two points are divisible by 5.
  // Those two are 5 * U(13*17).
  const std::vector<std::uint64_t> primes = {5, 5, 13, 17};
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
  const Quadruple mixed{imprimitive[0], imprimitive[1], primitive[0],
                        primitive[1]};
  EXPECT_TRUE(pte::cpu::is_primitive(points, primes, mixed));

  // A squarefree norm can never be imprimitive.
  const std::vector<std::uint64_t> sqfree = {5, 13, 17, 29};
  const auto sqfree_points = pte::cpu::point_set(sqfree);
  EXPECT_TRUE(pte::cpu::is_primitive(sqfree_points, sqfree, {0, 1, 2, 3}));
}

TEST(Enumerate, PrimitivityRejectsFullyImprimitiveQuadruples) {
  // Build a norm whose imprimitive family is large enough to hold all four
  // points: L1 = 5^2 * 13 * 17 * 29 has |U(13*17*29)| = 4 points scaled by 5.
  const std::vector<std::uint64_t> primes = {5, 5, 13, 17, 29};
  const auto points = pte::cpu::point_set(primes);

  std::vector<int> imprimitive;
  for (int i = 0; i < static_cast<int>(points.size()); ++i) {
    if (points[i].re % 5 == 0 && points[i].im % 5 == 0)
      imprimitive.push_back(i);
  }
  ASSERT_GE(imprimitive.size(), 4u);

  const Quadruple all_imprimitive{imprimitive[0], imprimitive[1],
                                  imprimitive[2], imprimitive[3]};
  EXPECT_FALSE(pte::cpu::is_primitive(points, primes, all_imprimitive));
}

TEST(Enumerate, ExemplarNormYieldsTheKnownSolution) {
  // The E_4 exemplar of paper.md section 7, reached by targeting its norm
  // rather than sweeping to 3.3e10.
  const std::vector<std::uint64_t> primes = {13, 89, 173, 233, 709};
  const auto points = pte::cpu::point_set(primes);
  ASSERT_EQ(points.size(), 16u);

  const auto nu = pte::cpu::nu_values(points);
  const auto hits = pte::cpu::meet_reference(nu);
  ASSERT_FALSE(hits.empty()) << "the known E_4 solution was not found";

  // The published solution has pair keys built from these nu values.
  const std::vector<std::int64_t> want = {11810488155, 32679910245, 15352988955,
                                          31172903205};
  bool matched = false;
  for (const Quadruple &q : hits) {
    std::vector<std::int64_t> got = {nu[q.a], nu[q.b], nu[q.c], nu[q.d]};
    std::vector<std::int64_t> expect = want;
    std::sort(got.begin(), got.end());
    std::sort(expect.begin(), expect.end());
    if (got == expect) {
      matched = true;
      EXPECT_TRUE(pte::cpu::verify_exact(nu, q));
      EXPECT_TRUE(pte::cpu::is_primitive(points, primes, q));
    }
  }
  EXPECT_TRUE(matched) << "hits found, but not the published quadruple";
}

} // namespace
