// Unit tests for the GMP-backed CPU big-integer code.
//
// The oracle is Boost.Multiprecision's cpp_int: an implementation with no
// shared code with GMP, so agreement is real evidence rather than a tautology.
#include "bigint.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <gtest/gtest.h>

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace {

using boost::multiprecision::cpp_int;
using pte::cpu::from_limbs;
using pte::cpu::kLimbBits;
using pte::cpu::limb_count;
using pte::cpu::to_limbs;

// The two libraries only need to agree on the value, so compare decimal
// strings.
std::string dec(const cpp_int &z) { return z.str(); }
std::string dec(const mpz_class &z) { return z.get_str(); }

// A handful of values that between them exercise carries, zero, ragged operand
// sizes and exact limb boundaries.
std::vector<std::string> InterestingDecimals() {
  std::vector<std::string> out = {
      "0",
      "1",
      "2",
      "4294967295",           // 2^32 - 1, one full limb of ones
      "4294967296",           // 2^32, carries into the next limb
      "18446744073709551615", // 2^64 - 1
      "123456789012345678901234567890",
  };
  // 2^k - 1 (all ones) and 2^k (single set bit) around several limb boundaries.
  for (int k : {31, 32, 33, 63, 64, 65, 127, 128, 255, 256}) {
    cpp_int p = cpp_int(1) << k;
    out.push_back(dec(p));
    out.push_back(dec(p - 1));
  }
  return out;
}

TEST(BigIntAdd, MatchesBoostCppInt) {
  const auto values = InterestingDecimals();
  for (const auto &a : values) {
    for (const auto &b : values) {
      const std::string expected = dec(cpp_int(a) + cpp_int(b));
      EXPECT_EQ(dec(pte::cpu::add(mpz_class(a), mpz_class(b))), expected)
          << "a=" << a << " b=" << b;
      EXPECT_EQ(pte::cpu::add_dec(a, b), expected) << "a=" << a << " b=" << b;
    }
  }
}

TEST(BigIntMul, MatchesBoostCppInt) {
  const auto values = InterestingDecimals();
  for (const auto &a : values) {
    for (const auto &b : values) {
      const std::string expected = dec(cpp_int(a) * cpp_int(b));
      EXPECT_EQ(dec(pte::cpu::mul(mpz_class(a), mpz_class(b))), expected)
          << "a=" << a << " b=" << b;
      EXPECT_EQ(pte::cpu::mul_dec(a, b), expected) << "a=" << a << " b=" << b;
    }
  }
}

TEST(BigIntMul, LopsidedOperands) {
  // 2^1000 - 1 times a single limb: every partial product carries.
  cpp_int big = (cpp_int(1) << 1000) - 1;
  const cpp_int small = 4294967295u;
  EXPECT_EQ(dec(pte::cpu::mul(mpz_class(dec(big)), mpz_class(dec(small)))),
            dec(big * small));
}

TEST(LimbCount, CountsLimbsNotBytes) {
  EXPECT_EQ(limb_count(mpz_class("0")), 0u);
  EXPECT_EQ(limb_count(mpz_class("1")), 1u);
  EXPECT_EQ(limb_count(mpz_class("4294967295")), 1u); // 2^32 - 1
  EXPECT_EQ(limb_count(mpz_class("4294967296")), 2u); // 2^32
  EXPECT_EQ(limb_count(mpz_class(dec((cpp_int(1) << 256) - 1))), 8u);
  EXPECT_EQ(limb_count(mpz_class(dec(cpp_int(1) << 256))), 9u);
}

TEST(Limbs, RoundTrips) {
  for (const auto &s : InterestingDecimals()) {
    const mpz_class z(s);
    // Round-trip at exactly the required width and at a padded width.
    for (std::size_t extra : {std::size_t{0}, std::size_t{1}, std::size_t{5}}) {
      const std::size_t nlimbs = limb_count(z) + extra;
      const auto limbs = to_limbs(z, nlimbs);
      ASSERT_EQ(limbs.size(), nlimbs);
      EXPECT_EQ(dec(from_limbs(limbs.data(), nlimbs)), s)
          << "value=" << s << " nlimbs=" << nlimbs;
    }
  }
}

TEST(Limbs, ZeroPadsEveryLimb) {
  // mpz_export writes nothing for zero; the padding loop has to cover the whole
  // buffer. Pre-fill with garbage so a missing fill would be caught.
  std::vector<std::uint32_t> limbs(8, 0xDEADBEEFu);
  to_limbs(mpz_class("0"), limbs.data(), limbs.size());
  for (std::size_t i = 0; i < limbs.size(); ++i)
    EXPECT_EQ(limbs[i], 0u) << "i=" << i;

  std::vector<std::uint32_t> partial(8, 0xDEADBEEFu);
  to_limbs(mpz_class("4294967296"), partial.data(), partial.size()); // 2^32
  EXPECT_EQ(partial[0], 0u);
  EXPECT_EQ(partial[1], 1u);
  for (std::size_t i = 2; i < partial.size(); ++i)
    EXPECT_EQ(partial[i], 0u) << "i=" << i;
}

TEST(Limbs, LittleEndianLimbOrder) {
  // 2^32 + 7 must land as limb0 = 7, limb1 = 1.
  const mpz_class z = (mpz_class(1) << kLimbBits) + 7;
  const auto limbs = to_limbs(z, 4);
  EXPECT_EQ(limbs[0], 7u);
  EXPECT_EQ(limbs[1], 1u);
  EXPECT_EQ(limbs[2], 0u);
  EXPECT_EQ(limbs[3], 0u);
}

TEST(Limbs, RejectsOverflowAndNegative) {
  const mpz_class two_pow_64 = mpz_class(1) << 64; // needs 3 limbs
  EXPECT_THROW(to_limbs(two_pow_64, 2), std::invalid_argument);
  EXPECT_NO_THROW(to_limbs(two_pow_64, 3));
  EXPECT_THROW(to_limbs(mpz_class("-1"), 4), std::invalid_argument);
}

} // namespace
