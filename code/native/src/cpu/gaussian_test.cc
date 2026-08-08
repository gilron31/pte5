// Unit tests for the Gaussian-integer layer.
#include "gaussian.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <set>
#include <vector>

namespace {

using pte::cpu::Gauss;

TEST(Gaussian, DecomposePrimeMatchesItsNorm) {
  // Every split prime below 200, plus a few larger ones. 1000000009 exercises
  // the 128-bit intermediates in the Tonelli-Shanks step.
  const std::vector<std::uint64_t> primes = {
      5,   13,  17,  29,  37,  41,      53,      61,      73,
      89,  97,  101, 109, 113, 137,     149,     157,     173,
      181, 193, 197, 233, 709, 1000033, 1000037, 1000081, 1000000009};
  for (const std::uint64_t p : primes) {
    Gauss z;
    ASSERT_TRUE(pte::cpu::decompose_prime(p, &z)) << "p = " << p;
    EXPECT_EQ(pte::cpu::norm(z), static_cast<std::int64_t>(p)) << "p = " << p;
    EXPECT_GE(z.re, z.im);
    EXPECT_GT(z.im, 0);
  }
}

TEST(Gaussian, DecomposePrimeRejectsNonSplit) {
  Gauss z;
  EXPECT_FALSE(pte::cpu::decompose_prime(2, &z));   // ramified
  EXPECT_FALSE(pte::cpu::decompose_prime(3, &z));   // inert
  EXPECT_FALSE(pte::cpu::decompose_prime(7, &z));   // inert
  EXPECT_FALSE(pte::cpu::decompose_prime(199, &z)); // inert
}

TEST(Gaussian, MulAndConjMulAgreeWithPlainMul) {
  const Gauss z{7, 3};
  const Gauss q{5, 2};
  Gauss a;
  Gauss b;
  pte::cpu::mul_and_conj_mul(z, q, &a, &b);
  EXPECT_EQ(a, pte::cpu::mul(z, q));
  EXPECT_EQ(b, pte::cpu::mul(z, pte::cpu::conj(q)));
}

TEST(Gaussian, CanonizeIsIdempotentAndUnitInvariant) {
  const Gauss z{3, 7};
  const Gauss c = pte::cpu::canonize(z);
  EXPECT_EQ(c, (Gauss{7, 3}));
  EXPECT_EQ(pte::cpu::canonize(c), c);
  // The whole orbit under units and conjugation collapses to one point.
  const std::vector<Gauss> orbit = {{3, 7},  {-7, 3}, {-3, -7}, {7, -3},
                                    {3, -7}, {7, 3},  {-3, 7},  {-7, -3}};
  for (const Gauss &w : orbit)
    EXPECT_EQ(pte::cpu::canonize(w), c);
}

TEST(Gaussian, PointSetHasTheExpectedSize) {
  // |U| = ceil(prod(b_i + 1) / 2).
  struct Case {
    std::vector<std::uint64_t> primes;
    std::size_t expected;
  };
  const std::vector<Case> cases = {
      {{5}, 1},             // 2/2
      {{5, 13}, 2},         // 4/2
      {{5, 13, 17}, 4},     // 8/2
      {{5, 5, 13, 17}, 6},  // 3*4/2
      {{5, 13, 17, 29}, 8}, // 16/2
      // The paper's exemplar norm, 13*89*173*233*709: five distinct primes.
      {{13, 89, 173, 233, 709}, 16},
  };
  for (const Case &c : cases) {
    const auto pts = pte::cpu::point_set(c.primes);
    EXPECT_EQ(pts.size(), c.expected);
    std::int64_t expected_norm = 1;
    for (const std::uint64_t p : c.primes)
      expected_norm *= static_cast<std::int64_t>(p);
    for (const Gauss &z : pts)
      EXPECT_EQ(pte::cpu::norm(z), expected_norm);
  }
}

TEST(Gaussian, NuIsInjectiveOnPointSets) {
  // The lemma the meet relies on: distinct canonized points of a common norm
  // never share nu, so a key collision cannot involve a repeated point.
  const std::vector<std::vector<std::uint64_t>> norms = {
      {5, 13, 17},
      {5, 13, 17, 29},
      {5, 5, 13, 17},
      {13, 17, 29, 37},
      {13, 89, 173, 233, 709}};
  for (const auto &primes : norms) {
    const auto pts = pte::cpu::point_set(primes);
    std::set<std::int64_t> seen;
    for (const Gauss &z : pts) {
      EXPECT_GE(pte::cpu::nu(z), 0);
      EXPECT_TRUE(seen.insert(pte::cpu::nu(z)).second) << "nu repeated";
    }
  }
}

TEST(Gaussian, ExemplarPointSetContainsTheKnownSolution) {
  // paper.md section 7: the E_4 exemplar. Its four root pairs come from four
  // points of norm L1 = 33065996717 via the (1+i) twist, e.g. the root pair
  // (A, B) = (252885, 46703) is the image of 149794 + 103091i.
  const auto pts = pte::cpu::point_set({13, 89, 173, 233, 709});
  ASSERT_EQ(pts.size(), 16u);

  const std::vector<Gauss> expected = {
      {149794, 103091}, {181309, 13894}, {155594, 94109}, {179219, 30766}};
  for (const Gauss &z : expected) {
    EXPECT_NE(std::find(pts.begin(), pts.end(), z), pts.end())
        << "missing point " << z.re << " + " << z.im << "i";
    // The twist really does produce the published root pair.
    EXPECT_EQ(pte::cpu::norm(z), 33065996717);
  }
  EXPECT_EQ(pte::cpu::nu(expected[0]), 11810488155); // = 252885 * 46703
}

} // namespace
