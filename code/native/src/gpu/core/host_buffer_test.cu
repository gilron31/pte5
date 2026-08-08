// Unit tests for the pinned host buffer.
//
// Correctness only. The reason HostBuffer exists -- pinned memory can be DMA'd
// directly, so it transfers faster and is what cudaMemcpyAsync needs to overlap
// at all -- is a bandwidth claim, and bandwidth is measured in
// apps/hybrid_demo.cc rather than asserted in a unit test.
#include "gpu/core/host_buffer.cuh"

#include "gpu/core/device.h"
#include "gpu/core/device_buffer.cuh"

#include <gtest/gtest.h>

#include <numeric>
#include <utility>
#include <vector>

namespace {

using pte::gpu::DeviceBuffer;
using pte::gpu::HostBuffer;

class GpuHostBuffer : public ::testing::Test {
protected:
  void SetUp() override {
    if (pte::gpu::device_count() == 0)
      GTEST_SKIP() << "no CUDA device available";
  }
};

TEST_F(GpuHostBuffer, RoundTripsThroughTheDevice) {
  constexpr std::size_t kCount = 4096;
  HostBuffer<int> src(kCount);
  ASSERT_EQ(src.count(), kCount);
  ASSERT_NE(src.data(), nullptr);

  // Fill the pinned memory directly, which is the point of exposing a range.
  std::iota(src.begin(), src.end(), 0);

  DeviceBuffer<int> dev(kCount);
  dev.CopyFromHost(src.data());

  HostBuffer<int> dst(kCount);
  dev.CopyToHost(dst.data());
  for (std::size_t i = 0; i < kCount; ++i)
    ASSERT_EQ(dst[i], static_cast<int>(i)) << "at " << i;
}

TEST_F(GpuHostBuffer, PartialCopyBackTouchesOnlyThePrefix) {
  constexpr std::size_t kCount = 256;
  HostBuffer<int> src(kCount);
  std::iota(src.begin(), src.end(), 1);

  DeviceBuffer<int> dev(kCount);
  dev.CopyFromHost(src.data());

  HostBuffer<int> dst(kCount);
  std::fill(dst.begin(), dst.end(), -1);
  dev.CopyToHost(dst.data(), 10);

  for (std::size_t i = 0; i < 10; ++i)
    EXPECT_EQ(dst[i], static_cast<int>(i) + 1);
  for (std::size_t i = 10; i < kCount; ++i)
    EXPECT_EQ(dst[i], -1) << "copied past the requested prefix";
}

TEST_F(GpuHostBuffer, EmptyBufferIsUsable) {
  HostBuffer<int> empty(0);
  EXPECT_EQ(empty.count(), 0u);
  EXPECT_EQ(empty.bytes(), 0u);
  EXPECT_EQ(empty.data(), nullptr);
  EXPECT_EQ(empty.begin(), empty.end());
}

TEST_F(GpuHostBuffer, IsMovable) {
  HostBuffer<int> a(16);
  std::iota(a.begin(), a.end(), 100);
  const int *raw = a.data();

  HostBuffer<int> b(std::move(a));
  EXPECT_EQ(b.data(), raw) << "move should transfer the allocation, not copy";
  EXPECT_EQ(b.count(), 16u);
  EXPECT_EQ(a.data(), nullptr); // NOLINT(bugprone-use-after-move)
  EXPECT_EQ(a.count(), 0u);     // NOLINT(bugprone-use-after-move)
  EXPECT_EQ(b[0], 100);

  HostBuffer<int> c(4);
  c = std::move(b);
  EXPECT_EQ(c.data(), raw);
  EXPECT_EQ(c.count(), 16u);
  EXPECT_EQ(c[15], 115);
}

} // namespace
