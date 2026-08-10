// The pool's contract is one allocation, not one per transfer -- that is the
// whole reason it exists (bench/gpu P7: 1.45 GB/s allocating per transfer
// against 13.4 GB/s pooled). So the load-bearing assertion here is on
// `allocations()`, not on any timing: a benchmark would tell us the machine is
// fast today, while this tells us the class is still doing the thing it was
// written for.
//
// A .cc, deliberately: cudaHostAlloc and cudaMemcpy are runtime-API calls, so
// this whole test compiles with g++ and gets -Wpedantic. If it ever needs nvcc,
// something device-side has leaked into a host header.
#include "gpu/core/pinned_pool.h"

#include "gpu/core/device.h"
#include "gpu/core/device_buffer.h"

#include <gtest/gtest.h>

#include <cstdint>
#include <numeric>
#include <vector>

namespace pte::gpu {
namespace {

class PinnedPoolTest : public ::testing::Test {
protected:
  void SetUp() override {
    if (device_count() == 0)
      GTEST_SKIP() << "no CUDA device";
  }
};

TEST_F(PinnedPoolTest, ReusesOneAllocationAcrossManyTransfers) {
  PinnedPool pool;
  DeviceBuffer<std::uint32_t> device(1024);
  std::vector<std::uint32_t> src(1024);
  std::iota(src.begin(), src.end(), 7u);

  pool.StageAndCopy(src.data(), src.size(), device);
  const std::size_t after_first = pool.allocations();
  EXPECT_EQ(after_first, 1u);

  for (int i = 0; i < 100; ++i)
    pool.StageAndCopy(src.data(), src.size(), device);
  EXPECT_EQ(pool.allocations(), after_first)
      << "the pool reallocated inside a steady-state loop, which is the P7 "
         "defect it exists to close";

  std::vector<std::uint32_t> back(1024);
  device.CopyToHost(back.data());
  EXPECT_EQ(back, src);
}

// Shrinking must not reallocate either: a per-window S-E pass hands over a
// different pair count every time, and a pool that reallocated on any change
// would be no better than cudaHostAlloc per call.
TEST_F(PinnedPoolTest, GrowsGeometricallyAndNeverShrinks) {
  PinnedPool pool;
  pool.Acquire<unsigned char>(4096);
  EXPECT_EQ(pool.allocations(), 1u);
  EXPECT_GE(pool.bytes(), 4096u);

  pool.Acquire<unsigned char>(1);
  pool.Acquire<unsigned char>(4096);
  EXPECT_EQ(pool.allocations(), 1u) << "shrinking must not reallocate";

  pool.Acquire<unsigned char>(4097);
  EXPECT_EQ(pool.allocations(), 2u);
  EXPECT_GE(pool.bytes(), 8192u) << "growth should be geometric, not exact";

  // A big jump is still one allocation, not a doubling per step.
  pool.Acquire<unsigned char>(10u << 20);
  EXPECT_EQ(pool.allocations(), 3u);
  EXPECT_GE(pool.bytes(), 10u << 20);
}

// Two ranges that must be live at once come from a single Acquire and a
// subspan; this is what the S-E window planner does with its lo/hi arrays, and
// it is the pattern the one-live-region contract is designed around.
TEST_F(PinnedPoolTest, OneAcquireCanBackTwoDisjointRanges) {
  PinnedPool pool;
  const std::span<std::uint32_t> all = pool.Acquire<std::uint32_t>(200);
  const std::span<std::uint32_t> lo = all.first(100);
  const std::span<std::uint32_t> hi = all.subspan(100);
  for (int i = 0; i < 100; ++i) {
    lo[static_cast<std::size_t>(i)] = static_cast<std::uint32_t>(i);
    hi[static_cast<std::size_t>(i)] = static_cast<std::uint32_t>(1000 + i);
  }
  EXPECT_EQ(lo[99], 99u);
  EXPECT_EQ(hi[0], 1000u);
  EXPECT_EQ(all[100], 1000u);
  EXPECT_EQ(pool.allocations(), 1u);
}

TEST_F(PinnedPoolTest, ZeroCountIsAnEmptySpanAndNoAllocation) {
  PinnedPool pool;
  EXPECT_TRUE(pool.Acquire<std::uint64_t>(0).empty());
  EXPECT_EQ(pool.allocations(), 0u);

  DeviceBuffer<std::uint32_t> device(4);
  pool.StageAndCopy<std::uint32_t>(nullptr, 0, device); // must not dereference
  EXPECT_EQ(pool.allocations(), 0u);
}

} // namespace
} // namespace pte::gpu
