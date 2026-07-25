// Unit tests for the GPU big-integer kernels.
//
// The oracle is a plain host schoolbook loop written here in the test, so this
// binary links pte_gpu and nothing else -- in particular no GMP. That is what
// makes the GPU suite independent of the CPU suite.
#include "gpu.h"

#include <gtest/gtest.h>

#include <cstdint>
#include <random>
#include <vector>

namespace {

using pte::gpu::kMaxLimbs;
using Limbs = std::vector<std::uint32_t>;

// Skips the whole binary when there is no usable device, so `ctest -L gpu`
// degrades to "skipped" rather than "failed" on a machine without a GPU.
class GpuBigInt : public ::testing::Test {
protected:
  void SetUp() override {
    if (pte::gpu::device_count() == 0)
      GTEST_SKIP() << "no CUDA device available";
  }
};

// Reference schoolbook multiply, host side, same limb layout.
Limbs HostMul(const Limbs &x, const Limbs &y) {
  const int n = static_cast<int>(x.size());
  Limbs acc(2 * n, 0);
  for (int i = 0; i < n; ++i) {
    std::uint32_t carry = 0;
    for (int j = 0; j < n; ++j) {
      const std::uint64_t p =
          static_cast<std::uint64_t>(x[i]) * y[j] + acc[i + j] + carry;
      acc[i + j] = static_cast<std::uint32_t>(p);
      carry = static_cast<std::uint32_t>(p >> 32);
    }
    acc[i + n] = carry;
  }
  return acc;
}

// Reference add with carry-out.
Limbs HostAdd(const Limbs &x, const Limbs &y) {
  const int n = static_cast<int>(x.size());
  Limbs out(n + 1, 0);
  std::uint32_t carry = 0;
  for (int i = 0; i < n; ++i) {
    const std::uint64_t s = static_cast<std::uint64_t>(x[i]) + y[i] + carry;
    out[i] = static_cast<std::uint32_t>(s);
    carry = static_cast<std::uint32_t>(s >> 32);
  }
  out[n] = carry;
  return out;
}

// The operand patterns worth covering at a given width.
std::vector<Limbs> Patterns(int nlimbs) {
  std::vector<Limbs> out;
  out.push_back(Limbs(nlimbs, 0));           // zero
  out.push_back(Limbs(nlimbs, 0xFFFFFFFFu)); // all ones: worst-case carries
  Limbs one(nlimbs, 0);
  one[0] = 1;
  out.push_back(one);
  Limbs high(nlimbs, 0);
  high[nlimbs - 1] = 0x80000000u; // single top bit
  out.push_back(high);
  Limbs low_ones(nlimbs, 0);
  low_ones[0] = 0xFFFFFFFFu; // one full limb of ones, rest zero
  out.push_back(low_ones);
  Limbs alt(nlimbs, 0);
  for (int i = 0; i < nlimbs; ++i)
    alt[i] = (i % 2 == 0) ? 0xAAAAAAAAu : 0x55555555u;
  out.push_back(alt);
  return out;
}

Limbs RandomLimbs(int nlimbs, std::mt19937 &rng) {
  std::uniform_int_distribution<std::uint32_t> dist(0, 0xFFFFFFFFu);
  Limbs v(nlimbs);
  for (auto &x : v)
    x = dist(rng);
  return v;
}

// Flattens a batch of operands into one contiguous buffer.
Limbs Flatten(const std::vector<Limbs> &rows) {
  Limbs out;
  for (const auto &r : rows)
    out.insert(out.end(), r.begin(), r.end());
  return out;
}

TEST_F(GpuBigInt, MulMatchesHostForAllPatternPairs) {
  for (int nlimbs : {1, 2, 4, 8}) {
    const auto pats = Patterns(nlimbs);
    // One batch entry per (a, b) pattern pair, all submitted together.
    std::vector<Limbs> as, bs;
    for (const auto &a : pats) {
      for (const auto &b : pats) {
        as.push_back(a);
        bs.push_back(b);
      }
    }
    const int batch = static_cast<int>(as.size());
    const Limbs flat_a = Flatten(as);
    const Limbs flat_b = Flatten(bs);
    Limbs got(static_cast<std::size_t>(batch) * 2 * nlimbs);

    pte::gpu::mul_batch(flat_a.data(), flat_b.data(), got.data(), nlimbs,
                        batch);

    for (int t = 0; t < batch; ++t) {
      const Limbs expected = HostMul(as[t], bs[t]);
      for (int i = 0; i < 2 * nlimbs; ++i) {
        ASSERT_EQ(got[static_cast<std::size_t>(t) * 2 * nlimbs + i],
                  expected[i])
            << "nlimbs=" << nlimbs << " batch item " << t << " limb " << i;
      }
    }
  }
}

TEST_F(GpuBigInt, MulMatchesHostOnRandomBatch) {
  std::mt19937 rng(1234);
  for (int nlimbs : {1, 3, 8, kMaxLimbs}) {
    const int batch = 1000; // deliberately not a multiple of the block size
    std::vector<Limbs> as, bs;
    for (int t = 0; t < batch; ++t) {
      as.push_back(RandomLimbs(nlimbs, rng));
      bs.push_back(RandomLimbs(nlimbs, rng));
    }
    const Limbs flat_a = Flatten(as);
    const Limbs flat_b = Flatten(bs);
    Limbs got(static_cast<std::size_t>(batch) * 2 * nlimbs);

    pte::gpu::mul_batch(flat_a.data(), flat_b.data(), got.data(), nlimbs,
                        batch);

    for (int t = 0; t < batch; ++t) {
      const Limbs expected = HostMul(as[t], bs[t]);
      for (int i = 0; i < 2 * nlimbs; ++i) {
        ASSERT_EQ(got[static_cast<std::size_t>(t) * 2 * nlimbs + i],
                  expected[i])
            << "nlimbs=" << nlimbs << " batch item " << t << " limb " << i;
      }
    }
  }
}

TEST_F(GpuBigInt, MulOfAllOnesIsExact) {
  // (2^64 - 1)^2 = 2^128 - 2^65 + 1, i.e. limbs 1,0,fffffffe,ffffffff.
  const Limbs a(2, 0xFFFFFFFFu);
  Limbs got(4);
  pte::gpu::mul_batch(a.data(), a.data(), got.data(), 2, 1);
  EXPECT_EQ(got[0], 0x00000001u);
  EXPECT_EQ(got[1], 0x00000000u);
  EXPECT_EQ(got[2], 0xFFFFFFFEu);
  EXPECT_EQ(got[3], 0xFFFFFFFFu);
}

TEST_F(GpuBigInt, AddMatchesHostForAllPatternPairs) {
  for (int nlimbs : {1, 2, 8, kMaxLimbs}) {
    const auto pats = Patterns(nlimbs);
    std::vector<Limbs> as, bs;
    for (const auto &a : pats) {
      for (const auto &b : pats) {
        as.push_back(a);
        bs.push_back(b);
      }
    }
    const int batch = static_cast<int>(as.size());
    const Limbs flat_a = Flatten(as);
    const Limbs flat_b = Flatten(bs);
    Limbs got(static_cast<std::size_t>(batch) * (nlimbs + 1));

    pte::gpu::add_batch(flat_a.data(), flat_b.data(), got.data(), nlimbs,
                        batch);

    for (int t = 0; t < batch; ++t) {
      const Limbs expected = HostAdd(as[t], bs[t]);
      for (int i = 0; i <= nlimbs; ++i) {
        ASSERT_EQ(got[static_cast<std::size_t>(t) * (nlimbs + 1) + i],
                  expected[i])
            << "nlimbs=" << nlimbs << " batch item " << t << " limb " << i;
      }
    }
  }
}

TEST_F(GpuBigInt, AddCarriesOutOfTheTopLimb) {
  const Limbs a(2, 0xFFFFFFFFu); // 2^64 - 1
  Limbs one(2, 0);
  one[0] = 1;
  Limbs got(3);
  pte::gpu::add_batch(a.data(), one.data(), got.data(), 2, 1);
  EXPECT_EQ(got[0], 0u);
  EXPECT_EQ(got[1], 0u);
  EXPECT_EQ(got[2], 1u); // the carry-out limb
}

TEST_F(GpuBigInt, EmptyBatchIsANoOp) {
  Limbs got(1, 0xDEADBEEFu);
  EXPECT_NO_THROW(pte::gpu::mul_batch(nullptr, nullptr, nullptr, 4, 0));
  EXPECT_NO_THROW(pte::gpu::add_batch(nullptr, nullptr, nullptr, 4, 0));
  EXPECT_EQ(got[0], 0xDEADBEEFu);
}

TEST(GpuBigIntArgs, RejectsOutOfRangeWidths) {
  // Argument validation happens before any device work, so this runs even with
  // no GPU present and is intentionally not a GpuBigInt test.
  Limbs dummy(kMaxLimbs, 0);
  Limbs out(4 * kMaxLimbs, 0);
  EXPECT_THROW(
      pte::gpu::mul_batch(dummy.data(), dummy.data(), out.data(), 0, 1),
      std::invalid_argument);
  EXPECT_THROW(
      pte::gpu::mul_batch(dummy.data(), dummy.data(), out.data(), -1, 1),
      std::invalid_argument);
  EXPECT_THROW(pte::gpu::mul_batch(dummy.data(), dummy.data(), out.data(),
                                   kMaxLimbs + 1, 1),
               std::invalid_argument);
  EXPECT_THROW(
      pte::gpu::add_batch(dummy.data(), dummy.data(), out.data(), 0, 1),
      std::invalid_argument);
  EXPECT_THROW(
      pte::gpu::mul_batch(dummy.data(), dummy.data(), out.data(), 4, -1),
      std::invalid_argument);
}

} // namespace
