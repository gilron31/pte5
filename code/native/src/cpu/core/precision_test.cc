// Tests for the precision policy.
//
// The point of this file is the wraparound fixture below. Everything else in
// the project agrees that the GPU meet's uint64 comparison is "a sound filter
// with false positives"; nothing until now demonstrated that a false positive
// exists, or that anything rejects one. A property no test can reach by
// accident -- a random collision mod 2^64 has probability 2^-64 -- is a
// property nobody has actually checked.
#include "cpu/core/precision.h"

#include <gtest/gtest.h>

#include <cstdint>

namespace {

using pte::cpu::exact_key_t;
using pte::cpu::kExactKeyBound;
using pte::cpu::kMaxNorm;
using pte::cpu::norm_t;

// Same opt-in the rest of the project uses; keeps -Wpedantic quiet about the
// GNU extension.
__extension__ using u128 = unsigned __int128;

// The same key the GPU kernel computes: modulo 2^64 by unsigned wraparound.
std::uint64_t TruncatedKey(std::uint64_t x, std::uint64_t y) {
  return x * x + y * y;
}

// The exact key, in the layer-2 arbiter type.
exact_key_t ExactKey(std::uint64_t x, std::uint64_t y) {
  return static_cast<exact_key_t>(x) * x + static_cast<exact_key_t>(y) * y;
}

TEST(Precision, ExactKeyBoundIsTheLastExactNorm) {
  // 2*L1^2 must fit in an unsigned 64-bit word at the bound...
  const auto k =
      static_cast<u128>(kExactKeyBound) * kExactKeyBound * 2;
  EXPECT_LE(k, static_cast<u128>(UINT64_MAX));

  // ...and must not, one above it.
  const auto over = static_cast<u128>(kExactKeyBound + 1) *
                    (kExactKeyBound + 1) * 2;
  EXPECT_GT(over, static_cast<u128>(UINT64_MAX));

  EXPECT_TRUE(pte::cpu::meet_is_exact(kExactKeyBound));
  EXPECT_FALSE(pte::cpu::meet_is_exact(kExactKeyBound + 1));
  EXPECT_FALSE(pte::cpu::meet_is_exact(kMaxNorm))
      << "the design bound is far above the exact regime, so confirmation is "
         "load-bearing for any real campaign";
}

// The fixture. Four distinct nu values whose two pair keys are different in
// exact arithmetic but congruent modulo 2^64 -- so the kernel emits this
// quadruple and the confirmation step is the only thing that rejects it.
//
//   nu       = [2^33, 3, 2^32, 2^63 - 3]
//   k(a,b)   = 2^66 + 9                      -> 9 mod 2^64
//   k(c,d)   = 2^64 + 2^126 - 6*2^63 + 9     -> 9 mod 2^64
TEST(Precision, WraparoundProducesAReachableFalsePositive) {
  const std::uint64_t a = std::uint64_t{1} << 33;
  const std::uint64_t b = 3;
  const std::uint64_t c = std::uint64_t{1} << 32;
  const std::uint64_t d = (std::uint64_t{1} << 63) - 3;

  // Four genuinely distinct points, so this is not a degenerate match.
  EXPECT_NE(a, c);
  EXPECT_NE(b, d);

  EXPECT_EQ(TruncatedKey(a, b), TruncatedKey(c, d))
      << "the kernel's comparison must accept this quadruple";
  EXPECT_EQ(TruncatedKey(a, b), 9u) << "both reduce to 9 mod 2^64";

  EXPECT_NE(ExactKey(a, b), ExactKey(c, d))
      << "but the exact keys differ, so it is a false positive";
}

// The direction that must never fail: a real solution is never lost. This is
// the half of the soundness argument that an exclusion claim depends on, since
// a false negative is invisible in the output by construction.
TEST(Precision, WraparoundHasNoFalseNegatives) {
  // Equal exact keys are always congruent, for any inputs whatever.
  for (std::uint64_t s = 1; s < 5000; s += 7) {
    // 1^2 + 8^2 == 4^2 + 7^2, scaled to sweep magnitudes including past 2^32
    // where the products themselves wrap.
    const std::uint64_t scale = s * 2654435761u; // spreads across the range
    const std::uint64_t p = 1 * scale, q = 8 * scale;
    const std::uint64_t r = 4 * scale, t = 7 * scale;
    ASSERT_EQ(ExactKey(p, q) == ExactKey(r, t), true) << "fixture broken";
    EXPECT_EQ(TruncatedKey(p, q), TruncatedKey(r, t))
        << "a true collision was lost by the filter at scale " << scale;
  }
}

TEST(Precision, ExactKeyTypeHoldsTheDesignBound) {
  // 2*B^2 at B = 10^12 is 81 bits; check it round-trips rather than trusting
  // the static_assert's sign test alone.
  const exact_key_t k = ExactKey(static_cast<std::uint64_t>(kMaxNorm),
                                 static_cast<std::uint64_t>(kMaxNorm));
  EXPECT_GT(k, 0);
  EXPECT_EQ(k / kMaxNorm / kMaxNorm, 2);
}

TEST(Precision, NormTypeHoldsTheDesignBound) {
  const norm_t n = kMaxNorm;
  EXPECT_EQ(n, 1'000'000'000'000);
}

} // namespace
