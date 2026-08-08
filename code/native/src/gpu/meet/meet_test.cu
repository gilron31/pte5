// Unit tests for the E_4 meet.
//
// The oracle is a host quadruple loop written here, so this binary depends on
// pte_gpu alone -- no GMP, no pte_cpu. That independence is the point: the CPU
// enumerator must never be able to hide a bug in the kernel.
//
// Tests go through the plain-C++ meet.h entry point, not the device-side
// signature, so they exercise the validation and staging as well as the kernel.
#include "gpu/meet/meet.h"

#include "gpu/core/device.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

namespace {

using pte::gpu::index_t;
using pte::gpu::MeetHit;
using pte::gpu::nu_t;

bool Less(const MeetHit &x, const MeetHit &y) {
  return std::tie(x.norm, x.a, x.b, x.c, x.d) <
         std::tie(y.norm, y.a, y.b, y.c, y.d);
}

class GpuMeet : public ::testing::Test {
protected:
  void SetUp() override {
    if (pte::gpu::device_count() == 0)
      GTEST_SKIP() << "no CUDA device available";
  }
};

// Same wraparound semantics as the kernel.
std::uint64_t PairKey(nu_t x, nu_t y) { return x * x + y * y; }

std::vector<MeetHit> HostMeet(const std::vector<nu_t> &nu,
                              const std::vector<index_t> &offsets) {
  std::vector<MeetHit> hits;
  for (index_t b = 0; b + 1 < static_cast<index_t>(offsets.size()); ++b) {
    const index_t begin = offsets[b];
    const int n = static_cast<int>(offsets[b + 1] - begin);
    const nu_t *v = nu.data() + begin;
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
  std::sort(hits.begin(), hits.end(), Less);
  return hits;
}

// Runs the meet and returns its hits, sorted so they can be compared with the
// oracle's regardless of the order the atomics produced.
std::vector<MeetHit> DeviceMeet(const std::vector<nu_t> &nu,
                                const std::vector<index_t> &offsets,
                                index_t max_hits, index_t *reported) {
  const index_t num_norms = static_cast<index_t>(offsets.size()) - 1;
  std::vector<MeetHit> hits(static_cast<std::size_t>(std::max<index_t>(max_hits, 1)));
  const index_t found = pte::gpu::meet_batch(
      nu.data(), static_cast<index_t>(nu.size()), offsets.data(), num_norms,
      hits.data(), max_hits);
  if (reported != nullptr)
    *reported = found;
  hits.resize(static_cast<std::size_t>(std::min(found, max_hits)));
  std::sort(hits.begin(), hits.end(), Less);
  return hits;
}

void ExpectMatchesHost(const std::vector<nu_t> &nu,
                       const std::vector<index_t> &offsets) {
  const auto expected = HostMeet(nu, offsets);
  index_t reported = 0;
  const auto got =
      DeviceMeet(nu, offsets, static_cast<index_t>(expected.size()) + 16,
                 &reported);
  EXPECT_EQ(reported, static_cast<index_t>(expected.size()));
  EXPECT_EQ(got, expected);
}

TEST_F(GpuMeet, FindsAPlantedCollision) {
  // 1^2 + 8^2 == 4^2 + 7^2.
  const std::vector<nu_t> nu = {1, 4, 7, 8, 11};
  const std::vector<index_t> offsets = {0, 5};
  index_t reported = 0;
  const auto got = DeviceMeet(nu, offsets, 8, &reported);
  ASSERT_EQ(reported, 1);
  ASSERT_EQ(got.size(), 1u);
  EXPECT_EQ(got[0], (MeetHit{0, 0, 3, 1, 2}));
  ExpectMatchesHost(nu, offsets);
}

TEST_F(GpuMeet, FindsNothingWhenThereIsNothing) {
  const std::vector<nu_t> nu = {1, 2, 4, 8, 16, 32};
  const std::vector<index_t> offsets = {0, 6};
  index_t reported = -1;
  EXPECT_TRUE(DeviceMeet(nu, offsets, 8, &reported).empty());
  EXPECT_EQ(reported, 0);
}

TEST_F(GpuMeet, HandlesEmptyAndUndersizedNorms) {
  EXPECT_EQ(pte::gpu::meet_batch(nullptr, 0, nullptr, 0, nullptr, 0), 0);

  // Norms with fewer than four points can never produce a hit.
  const std::vector<nu_t> nu = {3, 5, 7, 11, 13};
  const std::vector<index_t> offsets = {0, 0, 3, 5};
  index_t reported = -1;
  EXPECT_TRUE(DeviceMeet(nu, offsets, 8, &reported).empty());
  EXPECT_EQ(reported, 0);
}

TEST_F(GpuMeet, ManyNormsAtOnceMatchTheHost) {
  // A batch mixing sizes, with collisions planted in a few of them.
  std::mt19937_64 rng(20260808);
  std::vector<nu_t> nu;
  std::vector<index_t> offsets = {0};
  for (int b = 0; b < 500; ++b) {
    const int n = 4 + static_cast<int>(rng() % 9); // D between 4 and 12
    for (int t = 0; t < n; ++t)
      nu.push_back(rng() % 100000);
    if (b % 7 == 0) {
      // Force nu[0]^2 + nu[3]^2 == nu[1]^2 + nu[2]^2 via 1,4,7,8 scaled.
      const std::size_t base = nu.size() - n;
      const nu_t s = 1 + rng() % 1000;
      nu[base + 0] = 1 * s;
      nu[base + 1] = 4 * s;
      nu[base + 2] = 7 * s;
      nu[base + 3] = 8 * s;
    }
    offsets.push_back(static_cast<index_t>(nu.size()));
  }
  ExpectMatchesHost(nu, offsets);
}

TEST_F(GpuMeet, ReportsTruncationWithoutOverrunning) {
  // Four norms, each carrying one planted collision, but room for only two.
  std::vector<nu_t> nu;
  std::vector<index_t> offsets = {0};
  for (int b = 0; b < 4; ++b) {
    const nu_t s = static_cast<nu_t>(b) + 1;
    for (const nu_t base : {nu_t{1}, nu_t{4}, nu_t{7}, nu_t{8}})
      nu.push_back(base * s);
    offsets.push_back(static_cast<index_t>(nu.size()));
  }

  constexpr index_t kRoom = 2;
  std::vector<MeetHit> hits(kRoom, MeetHit{-1, -1, -1, -1, -1});
  const index_t found =
      pte::gpu::meet_batch(nu.data(), static_cast<index_t>(nu.size()),
                           offsets.data(), 4, hits.data(), kRoom);
  EXPECT_EQ(found, 4) << "the count must report every hit, not just the stored";
  // Exactly the first kRoom slots were filled, and with real hits.
  for (index_t t = 0; t < kRoom; ++t) {
    EXPECT_GE(hits[t].norm, 0);
    EXPECT_LT(hits[t].norm, 4);
  }
}

TEST_F(GpuMeet, RejectsBadArguments) {
  const std::vector<nu_t> nu = {1, 4, 7, 8};
  const std::vector<index_t> offsets = {0, 4};
  std::vector<MeetHit> hits(1);
  const index_t n = static_cast<index_t>(nu.size());

  EXPECT_THROW(
      pte::gpu::meet_batch(nu.data(), n, offsets.data(), -1, hits.data(), 1),
      std::invalid_argument);
  EXPECT_THROW(
      pte::gpu::meet_batch(nu.data(), n, offsets.data(), 1, hits.data(), -1),
      std::invalid_argument);
  EXPECT_THROW(
      pte::gpu::meet_batch(nullptr, n, offsets.data(), 1, hits.data(), 1),
      std::invalid_argument);
  EXPECT_THROW(
      pte::gpu::meet_batch(nu.data(), n, offsets.data(), 1, nullptr, 1),
      std::invalid_argument);

  // Offsets running past the nu array would be an out-of-range device read.
  const std::vector<index_t> too_long = {0, 99};
  EXPECT_THROW(
      pte::gpu::meet_batch(nu.data(), n, too_long.data(), 1, hits.data(), 1),
      std::invalid_argument);

  const std::vector<index_t> not_sorted = {0, 3, 2};
  EXPECT_THROW(
      pte::gpu::meet_batch(nu.data(), n, not_sorted.data(), 2, hits.data(), 1),
      std::invalid_argument);
}

TEST_F(GpuMeet, LargePointSetMatchesTheHost) {
  // D = 40 exercises the quadruple loop far past the common D = 4 case.
  std::mt19937_64 rng(99);
  std::vector<nu_t> nu;
  for (int t = 0; t < 40; ++t)
    nu.push_back(rng() % 1000000);
  nu[0] = 3;
  nu[1] = 12;
  nu[2] = 21;
  nu[3] = 24; // 3^2+24^2 == 12^2+21^2
  const std::vector<index_t> offsets = {0, 40};
  ExpectMatchesHost(nu, offsets);
}

} // namespace
