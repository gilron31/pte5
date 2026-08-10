// Champion point sets against the enumeration's own construction.
//
// The two are separate implementations at separate widths (see
// champion_points.h for why), so the fixture that matters is that they agree
// wherever both are valid -- every norm small enough for the int32 path. Above
// that only the invariants are checkable, and they are: D = 2^(k-1), every
// point has norm L1 exactly, and nu is injective.
#include "cpu/champion/champion_points.h"

#include "cpu/bigint/bigint.h"
#include "cpu/champion/champion.h"
#include "cpu/enumerate/enumerate.h"
#include "cpu/gaussian/gaussian.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <set>
#include <vector>

namespace pte::cpu {
namespace {

TEST(ChampionPoints, AgreesWithTheInt32PointSetWhereBothAreValid) {
  const std::vector<std::vector<factor_t>> cases{
      {5, 13},
      {5, 13, 17},
      {13, 89, 173, 233, 709}, // paper.md section 7's exemplar norm
      {5, 13, 17, 29, 37, 41},
  };
  for (const std::vector<factor_t> &primes : cases) {
    const ChampionNorm champ = champion_norm(primes);

    std::vector<factor_t> as_u64(primes.begin(), primes.end());
    const std::vector<norm_t> want = nu_values(point_set(as_u64));

    std::vector<core::U128> want_u128;
    want_u128.reserve(want.size());
    for (const norm_t v : want)
      want_u128.push_back(core::u128_from(static_cast<std::uint64_t>(v)));
    std::sort(want_u128.begin(), want_u128.end(),
              [](const core::U128 &a, const core::U128 &b) {
                return core::cmp(a, b) < 0;
              });
    std::vector<core::U128> got = champ.nu;
    std::sort(got.begin(), got.end(),
              [](const core::U128 &a, const core::U128 &b) {
                return core::cmp(a, b) < 0;
              });

    ASSERT_EQ(got.size(), want_u128.size()) << "point counts differ";
    for (std::size_t i = 0; i < got.size(); ++i)
      EXPECT_TRUE(got[i] == want_u128[i]) << "nu differs at " << i;
  }
}

TEST(ChampionPoints, DIsTwoToTheKMinusOne) {
  for (int k = 2; k <= 12; ++k) {
    const ChampionNorm champ = champion_norm(smallest_split_primes(k));
    EXPECT_EQ(champ.nu.size(), std::size_t{1} << (k - 1))
        << "k = " << k << " distinct primes must give D = 2^(k-1)";
    // nu injective: the precondition champion_meet checks for.
    const std::set<core::U128> distinct(champ.nu.begin(), champ.nu.end());
    EXPECT_EQ(distinct.size(), champ.nu.size()) << "nu is not injective";
  }
}

// The whole reason this file exists: at 20 primes the coordinates are 61 bits
// and nu is 122, so nothing on the int32 path could have produced this.
TEST(ChampionPoints, ReachesChampionWidths) {
  const ChampionNorm champ = champion_norm(smallest_split_primes(20));
  EXPECT_EQ(champ.nu.size(), 1u << 19) << "20 primes must give D = 2^19";

  int max_bits = 0;
  for (const core::U128 &v : champ.nu)
    max_bits = std::max(max_bits, core::bit_length(v));
  EXPECT_GT(max_bits, 118) << "nu should reach ~122 bits at this L1";
  EXPECT_LE(max_bits, 128);

  // L1 itself, and the key width the meet will need.
  EXPECT_GT(mpz_sizeinbase(champ.l1.get_mpz_t(), 2), 118u);
  EXPECT_LT(mpz_sizeinbase(champ.l1.get_mpz_t(), 2), 123u);

  // Every point really has norm L1: m^2 + n^2 == L1 for the (m,n) that nu and
  // L1 determine. Checking the identity nu <= L1 and L1 - nu even is the cheap
  // half; the full recovery is done for a sample.
  const mpz_class l1 = champ.l1;
  for (std::size_t i = 0; i < champ.nu.size(); i += 4096) {
    const mpz_class nu = from_uint<4>(champ.nu[i]);
    // m^2 = (L1 + nu)/2, n^2 = (L1 - nu)/2 -- both must be perfect squares.
    ASSERT_TRUE((l1 + nu) % 2 == 0);
    const mpz_class m2 = (l1 + nu) / 2, n2 = (l1 - nu) / 2;
    mpz_class m, n;
    mpz_sqrt(m.get_mpz_t(), m2.get_mpz_t());
    mpz_sqrt(n.get_mpz_t(), n2.get_mpz_t());
    EXPECT_EQ(m * m, m2) << "point " << i << " is not on the norm circle";
    EXPECT_EQ(n * n, n2) << "point " << i << " is not on the norm circle";
    EXPECT_EQ(m * m + n * n, l1);
  }
}

// A champion point set feeding straight into the meet is the end-to-end path
// S-E runs; at 10 primes it is small enough for a test and already past the
// widths the int32 enumeration can represent (L1 ~ 1.6e15 > 10^12).
TEST(ChampionPoints, FeedsTheChampionMeetDirectly) {
  const ChampionNorm champ = champion_norm(smallest_split_primes(10));
  ASSERT_EQ(champ.nu.size(), 512u);
  ChampionStats stats;
  const auto hits =
      champion_meet(std::span<const core::U128>(champ.nu), {}, &stats);
  EXPECT_EQ(stats.pairs, 512ull * 511 / 2);
  // Whether this particular norm has a collision is a fact about it, not
  // something to assert; what must hold is that anything reported is real.
  for (const core::PairCollision &h : hits) {
    EXPECT_TRUE(
        core::add(core::square(champ.nu[h.a]), core::square(champ.nu[h.b])) ==
        core::add(core::square(champ.nu[h.c]), core::square(champ.nu[h.d])));
  }
}

TEST(ChampionPoints, RejectsRepeatedAndNonSplitPrimes) {
  EXPECT_THROW(champion_norm({5, 13, 5}), std::invalid_argument);
  EXPECT_THROW(champion_norm({5, 7}), NotASplitPrime); // 7 = 3 mod 4, inert
  EXPECT_THROW(champion_norm(smallest_split_primes(22)), std::invalid_argument);
}

TEST(ChampionPoints, SmallestSplitPrimesAreTheExpectedOnes) {
  const std::vector<factor_t> want{5,  13, 17, 29, 37,  41,  53,
                                   61, 73, 89, 97, 101, 109, 113};
  EXPECT_EQ(smallest_split_primes(14), want);
  EXPECT_TRUE(smallest_split_primes(0).empty());
}

} // namespace
} // namespace pte::cpu
