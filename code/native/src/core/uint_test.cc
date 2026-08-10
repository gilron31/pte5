// UInt against an independent big-integer implementation.
//
// The oracle is Boost.Multiprecision's fixed-width uint256_t, chosen for the
// same reason cpu/bigint's test uses cpp_int: it is a second implementation
// that shares no code with the one under test. Its unchecked fixed-precision
// arithmetic wraps modulo 2^256, which is exactly the contract uint.h states
// for add and sub, so the two are comparable on overflow as well as below it.
//
// What is actually at stake here is not big-integer arithmetic in general -- it
// is that S-E's window boundaries are decided by `cmp` and its keys by `add`.
// A comparison that is wrong on one input in 2^64 relocates a pair into the
// wrong window, where its partner is never seen. So the fixtures below are
// deliberately weighted towards carries, borrows and limb boundaries rather
// than towards random middles.
//
// This covers the *host* definition of every operation. The device definition
// is a different instruction sequence behind the same signature (inline PTX
// carry chains), and gpu/bigint/uint_device_test.cu checks the two agree --
// which is the obligation core/hostdev.h attaches to sharing a type across the
// seam at all.
#include "core/uint.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <gtest/gtest.h>

#include <cstdint>
#include <random>
#include <vector>

namespace pte::core {
namespace {

using Big = boost::multiprecision::uint256_t;

template <int N> Big to_big(const UInt<N> &a) {
  Big r = 0;
  for (int i = N - 1; i >= 0; --i) {
    r <<= 32;
    r |= Big(a.w[i]);
  }
  return r;
}

U256 from_big(Big v) {
  U256 r{};
  for (int i = 0; i < 8; ++i) {
    r.w[i] = static_cast<std::uint32_t>(v & Big(0xffffffffu));
    v >>= 32;
  }
  return r;
}

// Readability helper: a U256 written as four 64-bit chunks, low first.
U256 make(std::uint64_t a, std::uint64_t b, std::uint64_t c, std::uint64_t d) {
  const std::uint64_t v[4] = {a, b, c, d};
  U256 r{};
  for (int i = 0; i < 4; ++i) {
    r.w[2 * i] = static_cast<std::uint32_t>(v[i]);
    r.w[2 * i + 1] = static_cast<std::uint32_t>(v[i] >> 32);
  }
  return r;
}

// The interesting inputs: limb boundaries, all-ones patterns, and values one
// step either side of a carry.
std::vector<U256> Corpus() {
  const std::uint64_t m = 0xffffffffffffffffull;
  std::vector<U256> v{
      make(0, 0, 0, 0), make(1, 0, 0, 0), make(m, 0, 0, 0), make(0, 1, 0, 0),
      make(m, m, 0, 0), make(0, 0, 1, 0), make(m, m, m, 0), make(0, 0, 0, 1),
      make(m, m, m, m), make(m, m, m, m - 1), make(1, 0, 0, 1ull << 63),
      // 32-bit limb boundaries specifically, since the limbs are 32 bits and
      // the 64-bit patterns above straddle two of them at a time.
      make(0xffffffffull, 0, 0, 0), make(0x100000000ull, 0, 0, 0)};
  std::mt19937_64 rng(20260810);
  for (int k = 0; k < 64; ++k) {
    U256 x{};
    // Vary how many limbs are populated, so narrow values are well covered.
    const int limbs = 1 + (k % 8);
    for (int i = 0; i < limbs; ++i)
      x.w[i] = static_cast<std::uint32_t>(rng());
    v.push_back(x);
  }
  return v;
}

TEST(U256, RoundTripsThroughTheOracle) {
  for (const U256 &a : Corpus())
    EXPECT_EQ(from_big(to_big(a)), a);
}

TEST(U256, CompareMatchesTheOracle) {
  const auto corpus = Corpus();
  for (const U256 &a : corpus) {
    for (const U256 &b : corpus) {
      const Big x = to_big(a), y = to_big(b);
      const int want = x < y ? -1 : (x == y ? 0 : 1);
      EXPECT_EQ(cmp(a, b), want) << "cmp disagrees on " << x << " vs " << y;
      EXPECT_EQ(a < b, x < y);
      EXPECT_EQ(a == b, x == y);
    }
  }
}

// The failure this guards against is the classic one: comparing limb arrays
// from the least significant end. Such a cmp gets the great majority of random
// inputs right and is wrong exactly when the high limbs differ.
TEST(U256, CompareLooksAtTheHighLimbFirst) {
  const std::uint64_t m = 0xffffffffffffffffull;
  EXPECT_GT(cmp(make(m, 0, 0, 1), make(0, 0, 0, 0)), 0);
  EXPECT_GT(cmp(make(0, 0, 0, 1), make(m, m, m, 0)), 0); // 2^192 > 2^192 - 1
}

TEST(U256, AddMatchesTheOracle) {
  const auto corpus = Corpus();
  for (const U256 &a : corpus)
    for (const U256 &b : corpus)
      EXPECT_EQ(to_big(add(a, b)), Big(to_big(a) + to_big(b)))
          << "add disagrees";
}

TEST(U256, SubMatchesTheOracle) {
  const auto corpus = Corpus();
  for (const U256 &a : corpus)
    for (const U256 &b : corpus)
      EXPECT_EQ(to_big(sub(a, b)), Big(to_big(a) - to_big(b)))
          << "sub disagrees";
}

TEST(U256, AddAndSubAreInverse) {
  const auto corpus = Corpus();
  for (const U256 &a : corpus)
    for (const U256 &b : corpus)
      EXPECT_EQ(sub(add(a, b), b), a);
}

TEST(U256, ShiftsMatchTheOracle) {
  const auto corpus = Corpus();
  for (const U256 &a : corpus) {
    for (int n : {0, 1, 7, 63, 64, 65, 127, 128, 191, 192, 200, 255}) {
      EXPECT_EQ(to_big(shl(a, n)), Big(to_big(a) << n)) << "shl by " << n;
      EXPECT_EQ(to_big(shr(a, n)), Big(to_big(a) >> n)) << "shr by " << n;
    }
    // Out of range is defined here rather than left as the UB of a 64-bit
    // shift by 64, which is how a shift helper usually goes wrong.
    EXPECT_TRUE(is_zero(shl(a, 256)));
    EXPECT_TRUE(is_zero(shr(a, 256)));
    EXPECT_TRUE(is_zero(shl(a, -1)));
  }
}

TEST(U256, BitLengthMatchesTheOracle) {
  for (const U256 &a : Corpus()) {
    const Big x = to_big(a);
    int want = 0;
    for (Big t = x; t != 0; t >>= 1)
      ++want;
    EXPECT_EQ(bit_length(a), want);
  }
}

TEST(U256, MultiplyIsExactAt128Bits) {
  std::mt19937_64 rng(4242);
  const std::uint64_t m = 0xffffffffffffffffull;
  std::vector<U128> corpus{u128_from(0, 0),
                           u128_from(1, 0),
                           u128_from(m, 0),
                           u128_from(0, 1),
                           u128_from(m, m),
                           u128_from(m, m - 1),
                           u128_from(1ull << 63, 1ull << 63)};
  for (int k = 0; k < 64; ++k)
    corpus.push_back(u128_from(rng(), k % 2 ? rng() : 0));

  for (const U128 &a : corpus) {
    for (const U128 &b : corpus) {
      const Big want = to_big(a) * to_big(b);
      EXPECT_EQ(to_big(mul(a, b)), want) << "mul disagrees";
    }
    EXPECT_EQ(to_big(square(a)), Big(to_big(a) * to_big(a)));
  }
}

// The champion regime this type exists for: nu is ~122 bits, the key is ~244,
// and the whole point is that it does not overflow. 2 * (2^122)^2 = 2^245.
TEST(U256, ChampionKeysFitWithHeadroom) {
  const U128 nu = u128_from(0xffffffffffffffffull, (1ull << 58) - 1); // < 2^122
  const U256 sq = square(nu);
  EXPECT_LT(bit_length(sq), 245);
  const U256 key = add(sq, sq);
  EXPECT_LT(bit_length(key), 246);
  EXPECT_EQ(to_big(key), Big(to_big(sq) * 2)); // no wrap
}

// nu values in the exhaustive regime are 64-bit; the same type has to carry
// them without ceremony, because the small-D cross-check against cpu/meet runs
// through exactly this path.
TEST(U256, NarrowValuesAgreeWithNativeArithmetic) {
  std::mt19937_64 rng(7);
  for (int k = 0; k < 1000; ++k) {
    // 31-bit nu, so that the native x*x + y*y the oracle uses cannot itself
    // overflow -- the check is that U256 agrees with uint64 where uint64 is
    // still right, which is the regime cpu/core/precision.h covers.
    const std::uint64_t x = rng() >> 33, y = rng() >> 33;
    const U256 key = add(square(u128_from(x)), square(u128_from(y)));
    const std::uint64_t want = x * x + y * y;
    EXPECT_EQ(chunk64(key, 0), want);
    EXPECT_EQ(chunk64(key, 1), 0u);
    EXPECT_EQ(chunk64(key, 2), 0u);
    EXPECT_EQ(chunk64(key, 3), 0u);
  }
}

} // namespace
} // namespace pte::core
