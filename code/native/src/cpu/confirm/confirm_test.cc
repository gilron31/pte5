// Tests for tower recovery and the certificate.
//
// The anchor is ExemplarReproducesTheReviewedTower: the three tower values in
// paper.md's Appendix were produced by drafts/verify_claims.py (sympy), and
// this file reproduces them from the same eight roots through an entirely
// different implementation. The Appendix's TODO asks for exactly that check
// plus the identity and the ladder before the exemplar is cited.
#include "cpu/confirm/confirm.h"

#include "cpu/enumerate/enumerate.h"
#include "cpu/generate/generate.h"
#include "cpu/meet/meet.h"

#include <gtest/gtest.h>

#include <array>
#include <stdexcept>
#include <vector>

namespace {

using pte::cpu::Certificate;
using pte::cpu::E4Sol;
using pte::cpu::Gauss;
using pte::cpu::Ladder;
using pte::cpu::norm_t;
using pte::cpu::Tower;

// The exemplar of paper.md 7, as it appears in drafts/verify_claims.py: eight
// roots in four pairs, each pair summing in squares to 2*L1.
constexpr long kExemplarRoots[8] = {252885, 46703,  195203, 167415,
                                    249703, 61485,  209985, 148453};

// Recover the four Gaussian points from the root pairs: r = m - n, r' = m + n.
std::array<Gauss, 4> ExemplarPoints() {
  std::array<Gauss, 4> pts{};
  for (int k = 0; k < 4; ++k) {
    const long a = kExemplarRoots[2 * k];
    const long b = kExemplarRoots[2 * k + 1];
    const long hi = a > b ? a : b;
    const long lo = a > b ? b : a;
    pts[static_cast<std::size_t>(k)] = Gauss{(hi + lo) / 2, (hi - lo) / 2};
  }
  return pts;
}

TEST(Confirm, ExemplarReproducesTheReviewedTower) {
  const std::array<Gauss, 4> pts = ExemplarPoints();
  const Tower t = pte::cpu::recover_tower(pts);

  EXPECT_EQ(t.l1, 33'065'996'717LL);
  // Above kExactKeyBound, so this exemplar is on the filter path: the
  // enumeration would have compared its keys modulo 2^64.
  EXPECT_FALSE(pte::cpu::meet_is_exact(t.l1));

  EXPECT_EQ(t.l2, mpz_class("489628056848329146064"));
  EXPECT_EQ(t.l3,
            mpz_class("175480010455650701584492675662518592000000"));
  EXPECT_EQ(t.sqrt_l4,
            mpz_class("40042900368028062136207226327668992000000"));

  // The eight roots must come back out, since they are what the identity is
  // stated over.
  std::vector<long> got;
  for (const mpz_class &r : t.roots)
    got.push_back(r.get_si());
  std::sort(got.begin(), got.end());
  std::vector<long> want(std::begin(kExemplarRoots), std::end(kExemplarRoots));
  std::sort(want.begin(), want.end());
  EXPECT_EQ(got, want);
}

TEST(Confirm, ExemplarSatisfiesTheDegree16IdentityAndTheLadder) {
  const Certificate c = pte::cpu::certify(ExemplarPoints());
  EXPECT_TRUE(c.identity_holds)
      << "prod(x^2 - r_k^2) != (((x^2-L1)^2-L2)^2-L3)^2-L4";
  EXPECT_TRUE(c.ladder.sqrt_l4_divisible) << "C'_8/2 = 151351200 | sqrt(L4)";
  EXPECT_TRUE(c.ladder.plus_is_square);
  EXPECT_TRUE(c.ladder.minus_is_square);
  EXPECT_TRUE(c.ladder.plus_divisible);
  EXPECT_TRUE(c.ladder.minus_divisible);
  EXPECT_TRUE(c.valid()) << c.failure();
}

// The constants are quoted from [CMSV, Table 1] via paper.md 4, so check them
// against their factorizations rather than against themselves.
TEST(Confirm, LadderConstantsMatchTheirFactorizations) {
  const long c8 = 64L * 27 * 25 * 49 * 11 * 13; // 2^6 3^3 5^2 7^2 11 13
  EXPECT_EQ(c8, 302'702'400L);
  EXPECT_EQ(pte::cpu::kC8Half, c8 / 2);
  EXPECT_EQ(pte::cpu::kC4Half, (4L * 9) / 2); // C'_4 = 2^2 3^2 = 36
}

// nu^2 + u^2 = L1^2, so the meet's nu condition and the tower's u condition are
// the same equation. This is why L2 exists for every hit, and it is the one
// algebraic fact linking the search to the object being searched for.
TEST(Confirm, TheMeetConditionIsTheL2Condition) {
  const pte::cpu::Generator g(200'000);
  pte::cpu::MeetScratch<pte::cpu::trunc_key_t> scratch;
  std::size_t checked = 0;

  g.for_each_norm([&](const pte::cpu::NormView &norm) {
    for (const Gauss &z : norm.points) {
      const auto nu = static_cast<norm_t>(z.nu());
      const auto u = static_cast<norm_t>(2 * z.re * z.im);
      ASSERT_EQ(static_cast<pte::cpu::exact_key_t>(nu) * nu +
                    static_cast<pte::cpu::exact_key_t>(u) * u,
                static_cast<pte::cpu::exact_key_t>(norm.l1) * norm.l1)
          << "L1 = " << norm.l1;
    }

    pte::cpu::meet<pte::cpu::trunc_key_t>(
        norm.D(), norm.nu.data(), scratch, [&](const E4Sol &s) {
          // Every meet hit must peel: recover_tower throws if the u-form
          // condition fails, which would mean the two forms are not equivalent.
          const Certificate c =
              pte::cpu::certify(pte::cpu::quadruple(norm.points, s));
          EXPECT_TRUE(c.valid()) << c.failure();
          ++checked;
        });
  });
  EXPECT_GT(checked, 0u) << "a run that certified nothing would pass vacuously";
}

// Every solution the enumeration finds below 10^6 must carry a certificate.
// This is the property the campaign will assert on every catalog entry.
TEST(Confirm, EverySolutionToOneMillionCertifies) {
  const pte::cpu::Generator g(1'000'000);
  pte::cpu::MeetScratch<pte::cpu::trunc_key_t> scratch;
  std::vector<norm_t> found;

  g.for_each_norm([&](const pte::cpu::NormView &norm) {
    pte::cpu::meet<pte::cpu::trunc_key_t>(
        norm.D(), norm.nu.data(), scratch, [&](const E4Sol &s) {
          const Certificate c =
              pte::cpu::certify(pte::cpu::quadruple(norm.points, s));
          EXPECT_TRUE(c.valid()) << c.failure();
          EXPECT_GT(c.tower.l2, 0);
          found.push_back(norm.l1);
        });
  });
  EXPECT_EQ(found.size(), 2u) << "the known E_4 count below 10^6";
}

// A truncated-key false positive must not reach a certificate. recover_tower
// recomputes the meet condition in GMP, which is the layer-2 arbiter for
// anything the 64-bit filter emits.
TEST(Confirm, RejectsAQuadrupleThatDoesNotMeet) {
  // Four genuine points of norm 5525 = 5^2 * 13 * 17, in a pairing that does
  // not satisfy the condition.
  const std::vector<Gauss> pts = pte::cpu::point_set({5, 5, 13, 17});
  ASSERT_GE(pts.size(), 4u);
  const std::array<Gauss, 4> quad{pts[0], pts[1], pts[2], pts[3]};

  const auto key = [](const Gauss &x, const Gauss &y) {
    return static_cast<pte::cpu::exact_key_t>(x.nu()) * x.nu() +
           static_cast<pte::cpu::exact_key_t>(y.nu()) * y.nu();
  };
  ASSERT_NE(key(quad[0], quad[1]), key(quad[2], quad[3]))
      << "this fixture is only meaningful on a non-solution";

  EXPECT_THROW(pte::cpu::recover_tower(quad), std::invalid_argument);
}

TEST(Confirm, RejectsMalformedInput) {
  const std::array<Gauss, 4> mixed_norms{Gauss{2, 1}, Gauss{2, 1}, Gauss{3, 2},
                                         Gauss{2, 1}};
  EXPECT_THROW(pte::cpu::recover_tower(mixed_norms), std::invalid_argument);

  // Not canonized: m < n would make the first root of that pair negative.
  const std::array<Gauss, 4> uncanonized{Gauss{1, 2}, Gauss{2, 1}, Gauss{2, 1},
                                         Gauss{2, 1}};
  EXPECT_THROW(pte::cpu::recover_tower(uncanonized), std::invalid_argument);

  const std::vector<Gauss> pts = pte::cpu::point_set({5, 13, 17});
  EXPECT_THROW(pte::cpu::quadruple(pts, E4Sol{0, 1, 2, 99}),
               std::invalid_argument);
  EXPECT_THROW(pte::cpu::quadruple(pts, E4Sol{0, 1, 2, -1}),
               std::invalid_argument);
}

// The identity is checked coefficient by coefficient, so perturbing any single
// tower value must break it. Without this, "identity_holds" could be a function
// that returns true.
TEST(Confirm, IdentityIsSensitiveToEveryTowerValue) {
  const Tower good = pte::cpu::recover_tower(ExemplarPoints());
  ASSERT_TRUE(pte::cpu::verify_identity(good));

  Tower t = good;
  t.l2 += 1;
  EXPECT_FALSE(pte::cpu::verify_identity(t));

  t = good;
  t.l3 += 1;
  EXPECT_FALSE(pte::cpu::verify_identity(t));

  t = good;
  t.sqrt_l4 += 1;
  EXPECT_FALSE(pte::cpu::verify_identity(t));

  t = good;
  t.l1 += 1;
  EXPECT_FALSE(pte::cpu::verify_identity(t));

  t = good;
  t.roots[3] += 1;
  EXPECT_FALSE(pte::cpu::verify_identity(t));
}

} // namespace
