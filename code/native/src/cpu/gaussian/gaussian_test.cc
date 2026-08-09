// Unit tests for the Gaussian-integer layer.
#include "cpu/gaussian/gaussian.h"

#include <gmpxx.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <random>
#include <set>
#include <vector>

namespace {

using pte::cpu::Gauss;
using u64 = std::uint64_t;

TEST(Gaussian, DecomposePrimeMatchesItsNorm) {
  // Every split prime below 200, plus a few larger ones. 1000000009 exercises
  // the 128-bit intermediates in the Tonelli-Shanks step.
  const std::vector<std::uint64_t> primes = {
      5,   13,  17,  29,  37,  41,      53,      61,      73,
      89,  97,  101, 109, 113, 137,     149,     157,     173,
      181, 193, 197, 233, 709, 1000033, 1000037, 1000081, 1000000009};
  for (const std::uint64_t p : primes) {
    const Gauss z = pte::cpu::decompose_prime(p);
    EXPECT_EQ(z.norm(), static_cast<std::int64_t>(p)) << "p = " << p;
    EXPECT_GE(z.re, z.im);
    EXPECT_GT(z.im, 0);
  }
}

// It must *throw*, not report failure through a return value. A silent
// "no decomposition" would propagate as "this norm has no points", which is
// indistinguishable from "this norm has no solutions" -- the one failure mode
// an exclusion search cannot tolerate.
TEST(Gaussian, DecomposePrimeThrowsOnNonSplitInput) {
  EXPECT_THROW(pte::cpu::decompose_prime(2), pte::cpu::NotASplitPrime);
  EXPECT_THROW(pte::cpu::decompose_prime(3), pte::cpu::NotASplitPrime);
  EXPECT_THROW(pte::cpu::decompose_prime(7), pte::cpu::NotASplitPrime);
  EXPECT_THROW(pte::cpu::decompose_prime(199), pte::cpu::NotASplitPrime);
  EXPECT_THROW(pte::cpu::decompose_prime(1), pte::cpu::NotASplitPrime);
  EXPECT_THROW(pte::cpu::decompose_prime(0), pte::cpu::NotASplitPrime);
}

// Composites congruent to 1 mod 4 are the dangerous case: 25 = 3^2 + 4^2 is a
// sum of two squares, so a decomposition exists and the old code would have
// returned it happily for something that is not prime.
TEST(Gaussian, DecomposePrimeThrowsOnCompositesCongruentTo1Mod4) {
  for (const std::uint64_t n : {u64{25}, u64{65}, u64{85}, u64{125}, u64{145},
                               u64{221}, u64{1105}, u64{9}, u64{21}, u64{33}})
    EXPECT_THROW(pte::cpu::decompose_prime(n), pte::cpu::NotASplitPrime)
        << "n = " << n;
}

TEST(Gaussian, MulAndConjMulAgreeWithPlainMul) {
  const Gauss z{7, 3};
  const Gauss q{5, 2};
  Gauss a;
  Gauss b;
  z.mul_and_conj_mul(q, &a, &b);
  EXPECT_EQ(a, z * q);
  EXPECT_EQ(b, z * q.conj());
}

TEST(Gaussian, CanonizeIsIdempotentAndUnitInvariant) {
  const Gauss z{3, 7};
  const Gauss c = z.canonize();
  EXPECT_EQ(c, (Gauss{7, 3}));
  EXPECT_EQ(c.canonize(), c);
  // The whole orbit under units and conjugation collapses to one point.
  const std::vector<Gauss> orbit = {{3, 7},  {-7, 3}, {-3, -7}, {7, -3},
                                    {3, -7}, {7, 3},  {-3, 7},  {-7, -3}};
  for (const Gauss &w : orbit)
    EXPECT_EQ(w.canonize(), c);
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
      EXPECT_EQ(z.norm(), expected_norm);
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
      EXPECT_GE(z.nu(), 0);
      EXPECT_TRUE(seen.insert(z.nu()).second) << "nu repeated";
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
    EXPECT_EQ(z.norm(), 33065996717);
  }
  EXPECT_EQ(expected[0].nu(), 11810488155); // = 252885 * 46703
}


// D14: don't trust a hand-rolled isqrt. GMP is already a dependency, so the
// wheel's authority is available without paying its per-call cost -- check
// against mpz_sqrt here rather than calling it in the descent.
TEST(Gaussian, IsqrtAgreesWithGmp) {
  const auto check = [](std::uint64_t n) {
    mpz_class m(std::to_string(n));
    mpz_class root;
    mpz_sqrt(root.get_mpz_t(), m.get_mpz_t());
    EXPECT_EQ(pte::cpu::isqrt(n), root.get_ui()) << "n = " << n;
  };

  // Boundaries first: perfect squares and their neighbours are where a
  // correction loop that is off by one shows up.
  for (std::uint64_t r = 0; r < 3000; ++r) {
    const std::uint64_t sq = r * r;
    check(sq);
    if (sq > 0)
      check(sq - 1);
    check(sq + 1);
  }

  // Past 2^53 the double seed cannot represent n exactly, so both correction
  // directions get exercised.
  check(0);
  check(1);
  check(UINT64_MAX);
  check((std::uint64_t{1} << 53) + 1);
  for (std::uint64_t r = (std::uint64_t{1} << 31) - 4;
       r < (std::uint64_t{1} << 31) + 4; ++r)
    check(r * r);

  std::mt19937_64 rng(20260809);
  for (int i = 0; i < 20000; ++i)
    check(rng());
}

TEST(Gaussian, IsPrimeAgreesWithGmp) {
  const auto check = [](std::uint64_t n) {
    mpz_class m(std::to_string(n));
    // probab_prime_p returns 0 (composite), 1 (probably), 2 (definitely).
    const bool gmp_says = mpz_probab_prime_p(m.get_mpz_t(), 40) != 0;
    EXPECT_EQ(pte::cpu::is_prime(n), gmp_says) << "n = " << n;
  };
  for (std::uint64_t n = 0; n < 5000; ++n)
    check(n);
  // Straddle the 3,215,031,751 witness-set boundary in both directions.
  for (std::uint64_t n = 3215031700ull; n < 3215031800ull; ++n)
    check(n);
  std::mt19937_64 rng(99);
  for (int i = 0; i < 3000; ++i)
    check(rng() % 1000000000000ull);
}

} // namespace