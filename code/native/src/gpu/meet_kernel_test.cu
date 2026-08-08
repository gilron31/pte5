// Unit tests for the E_4 meet kernel.
//
// The oracle is a host quadruple loop written here, so this binary depends on
// pte_gpu alone -- no GMP, no pte_cpu. That independence is the point: the CPU
// enumerator must never be able to hide a bug in the kernel.
#include "gpu.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <random>
#include <vector>

namespace {

using Hit = std::array<int, 5>;

class GpuMeet : public ::testing::Test {
protected:
  void SetUp() override {
    if (pte::gpu::device_count() == 0)
      GTEST_SKIP() << "no CUDA device available";
  }
};

// Same wraparound semantics as the kernel.
std::uint64_t PairKey(std::uint64_t x, std::uint64_t y) {
  return x * x + y * y;
}

std::vector<Hit> HostMeet(const std::vector<std::uint64_t> &nu,
                          const std::vector<std::int32_t> &offsets) {
  std::vector<Hit> hits;
  for (int b = 0; b + 1 < static_cast<int>(offsets.size()); ++b) {
    const int begin = offsets[b];
    const int n = offsets[b + 1] - begin;
    const std::uint64_t *v = nu.data() + begin;
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j)
        for (int k = j + 1; k < n; ++k)
          for (int l = k + 1; l < n; ++l) {
            if (PairKey(v[i], v[j]) == PairKey(v[k], v[l]))
              hits.push_back({b, i, j, k, l});
            if (PairKey(v[i], v[k]) == PairKey(v[j], v[l]))
              hits.push_back({b, i, k, j, l});
            if (PairKey(v[i], v[l]) == PairKey(v[j], v[k]))
              hits.push_back({b, i, l, j, k});
          }
  }
  std::sort(hits.begin(), hits.end());
  return hits;
}

// Runs the kernel and returns its hits, sorted so they can be compared with the
// oracle's regardless of the order the atomics produced.
std::vector<Hit> DeviceMeet(const std::vector<std::uint64_t> &nu,
                            const std::vector<std::int32_t> &offsets,
                            int max_hits, int *reported) {
  const int num_norms = static_cast<int>(offsets.size()) - 1;
  std::vector<std::int32_t> raw(
      static_cast<std::size_t>(std::max(max_hits, 1)) *
      pte::gpu::kMeetHitStride);
  const int found = pte::gpu::meet_batch(nu.data(), offsets.data(), num_norms,
                                         raw.data(), max_hits);
  if (reported != nullptr)
    *reported = found;

  std::vector<Hit> hits;
  const int written = std::min(found, max_hits);
  for (int t = 0; t < written; ++t) {
    const std::size_t base =
        static_cast<std::size_t>(t) * pte::gpu::kMeetHitStride;
    hits.push_back({raw[base], raw[base + 1], raw[base + 2], raw[base + 3],
                    raw[base + 4]});
  }
  std::sort(hits.begin(), hits.end());
  return hits;
}

void ExpectMatchesHost(const std::vector<std::uint64_t> &nu,
                       const std::vector<std::int32_t> &offsets) {
  const auto expected = HostMeet(nu, offsets);
  int reported = 0;
  const auto got = DeviceMeet(
      nu, offsets, static_cast<int>(expected.size()) + 16, &reported);
  EXPECT_EQ(reported, static_cast<int>(expected.size()));
  EXPECT_EQ(got, expected);
}

TEST_F(GpuMeet, FindsAPlantedCollision) {
  // 1^2 + 8^2 == 4^2 + 7^2.
  const std::vector<std::uint64_t> nu = {1, 4, 7, 8, 11};
  const std::vector<std::int32_t> offsets = {0, 5};
  int reported = 0;
  const auto got = DeviceMeet(nu, offsets, 8, &reported);
  ASSERT_EQ(reported, 1);
  ASSERT_EQ(got.size(), 1u);
  EXPECT_EQ(got[0], (Hit{0, 0, 3, 1, 2}));
  ExpectMatchesHost(nu, offsets);
}

TEST_F(GpuMeet, FindsNothingWhenThereIsNothing) {
  const std::vector<std::uint64_t> nu = {1, 2, 4, 8, 16, 32};
  const std::vector<std::int32_t> offsets = {0, 6};
  int reported = -1;
  EXPECT_TRUE(DeviceMeet(nu, offsets, 8, &reported).empty());
  EXPECT_EQ(reported, 0);
}

TEST_F(GpuMeet, HandlesEmptyAndUndersizedNorms) {
  EXPECT_EQ(pte::gpu::meet_batch(nullptr, nullptr, 0, nullptr, 0), 0);

  // Norms with fewer than four points can never produce a hit.
  const std::vector<std::uint64_t> nu = {3, 5, 7, 11, 13};
  const std::vector<std::int32_t> offsets = {0, 0, 3, 5};
  int reported = -1;
  EXPECT_TRUE(DeviceMeet(nu, offsets, 8, &reported).empty());
  EXPECT_EQ(reported, 0);
}

TEST_F(GpuMeet, ManyNormsAtOnceMatchTheHost) {
  // A batch mixing sizes, with collisions planted in a few of them.
  std::mt19937_64 rng(20260808);
  std::vector<std::uint64_t> nu;
  std::vector<std::int32_t> offsets = {0};
  for (int b = 0; b < 500; ++b) {
    const int n = 4 + static_cast<int>(rng() % 9); // D between 4 and 12
    for (int t = 0; t < n; ++t)
      nu.push_back(rng() % 100000);
    if (b % 7 == 0) {
      // Force nu[0]^2 + nu[3]^2 == nu[1]^2 + nu[2]^2 via 1,4,7,8 scaled.
      const std::size_t base = nu.size() - n;
      const std::uint64_t s = 1 + rng() % 1000;
      nu[base + 0] = 1 * s;
      nu[base + 1] = 4 * s;
      nu[base + 2] = 7 * s;
      nu[base + 3] = 8 * s;
    }
    offsets.push_back(static_cast<std::int32_t>(nu.size()));
  }
  ExpectMatchesHost(nu, offsets);
}

TEST_F(GpuMeet, ReportsTruncationWithoutOverrunning) {
  // Four norms, each carrying one planted collision, but room for only two.
  std::vector<std::uint64_t> nu;
  std::vector<std::int32_t> offsets = {0};
  for (int b = 0; b < 4; ++b) {
    const std::uint64_t s = b + 1;
    for (const std::uint64_t base : {1ull, 4ull, 7ull, 8ull})
      nu.push_back(base * s);
    offsets.push_back(static_cast<std::int32_t>(nu.size()));
  }

  constexpr int kRoom = 2;
  std::vector<std::int32_t> raw(
      static_cast<std::size_t>(kRoom) * pte::gpu::kMeetHitStride, -1);
  const int found =
      pte::gpu::meet_batch(nu.data(), offsets.data(), 4, raw.data(), kRoom);
  EXPECT_EQ(found, 4) << "the count must report every hit, not just the stored";
  // Exactly the first kRoom slots were filled, and with real hits.
  for (int t = 0; t < kRoom; ++t) {
    const std::size_t base =
        static_cast<std::size_t>(t) * pte::gpu::kMeetHitStride;
    EXPECT_GE(raw[base], 0);
    EXPECT_LT(raw[base], 4);
  }
}

TEST_F(GpuMeet, RejectsBadArguments) {
  const std::vector<std::uint64_t> nu = {1, 4, 7, 8};
  const std::vector<std::int32_t> offsets = {0, 4};
  std::vector<std::int32_t> raw(pte::gpu::kMeetHitStride);
  EXPECT_THROW(
      pte::gpu::meet_batch(nu.data(), offsets.data(), -1, raw.data(), 1),
      std::invalid_argument);
  EXPECT_THROW(
      pte::gpu::meet_batch(nu.data(), offsets.data(), 1, raw.data(), -1),
      std::invalid_argument);
  EXPECT_THROW(pte::gpu::meet_batch(nullptr, offsets.data(), 1, raw.data(), 1),
               std::invalid_argument);
  EXPECT_THROW(pte::gpu::meet_batch(nu.data(), offsets.data(), 1, nullptr, 1),
               std::invalid_argument);
}

TEST_F(GpuMeet, LargePointSetMatchesTheHost) {
  // D = 40 exercises the quadruple loop far past the common D = 4 case.
  std::mt19937_64 rng(99);
  std::vector<std::uint64_t> nu;
  for (int t = 0; t < 40; ++t)
    nu.push_back(rng() % 1000000);
  nu[0] = 3;
  nu[1] = 12;
  nu[2] = 21;
  nu[3] = 24; // 3^2+24^2 == 12^2+21^2
  const std::vector<std::int32_t> offsets = {0, 40};
  ExpectMatchesHost(nu, offsets);
}

} // namespace
