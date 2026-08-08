// Unit tests for the append-only device output buffer.
//
// The overflow contract is the part worth pinning down: every append bumps the
// counter, including the ones that do not fit, so total() reports what the
// kernel wanted to write while only the first capacity() items were stored.
#include "gpu/core/device_append_buffer.cuh"

#include "gpu/core/device.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

namespace {

using pte::gpu::DeviceAppendBuffer;
using pte::gpu::DeviceAppendView;
using pte::gpu::index_t;

class GpuAppendBuffer : public ::testing::Test {
protected:
  void SetUp() override {
    if (pte::gpu::device_count() == 0)
      GTEST_SKIP() << "no CUDA device available";
  }
};

__global__ void PushRange(DeviceAppendView<int> out, int count) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < count)
    out.push(i);
}

std::vector<int> Drain(DeviceAppendBuffer<int> &buf, index_t *total) {
  const auto found = buf.total();
  if (total != nullptr)
    *total = static_cast<index_t>(found);
  const index_t stored = buf.stored(found);
  std::vector<int> out(static_cast<std::size_t>(stored));
  buf.CopyToHost(out.data(), stored);
  std::sort(out.begin(), out.end());
  return out;
}

TEST_F(GpuAppendBuffer, StoresEverythingWhenItFits) {
  DeviceAppendBuffer<int> buf(64);
  PushRange<<<1, 64>>>(buf.view(), 10);
  ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

  index_t total = 0;
  const auto got = Drain(buf, &total);
  EXPECT_EQ(total, 10);
  ASSERT_EQ(got.size(), 10u);
  for (int i = 0; i < 10; ++i)
    EXPECT_EQ(got[i], i);
}

TEST_F(GpuAppendBuffer, CountsPastCapacityButStoresOnlyCapacity) {
  constexpr index_t kCapacity = 8;
  DeviceAppendBuffer<int> buf(kCapacity);
  PushRange<<<1, 128>>>(buf.view(), 100);
  ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

  index_t total = 0;
  const auto got = Drain(buf, &total);
  EXPECT_EQ(total, 100) << "the counter must report every attempted append";
  EXPECT_EQ(static_cast<index_t>(got.size()), kCapacity);
  // Whatever landed is a real value pushed by some thread, not garbage.
  for (const int v : got) {
    EXPECT_GE(v, 0);
    EXPECT_LT(v, 100);
  }
  EXPECT_EQ(buf.capacity(), kCapacity);
}

TEST_F(GpuAppendBuffer, EmptyBufferIsUsable) {
  DeviceAppendBuffer<int> buf(0);
  PushRange<<<1, 32>>>(buf.view(), 5);
  ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

  index_t total = 0;
  EXPECT_TRUE(Drain(buf, &total).empty());
  EXPECT_EQ(total, 5) << "counting must work even with nowhere to store";
}

TEST_F(GpuAppendBuffer, StartsEmpty) {
  DeviceAppendBuffer<int> buf(16);
  EXPECT_EQ(buf.total(), 0u);
  EXPECT_EQ(buf.stored(0), 0);
}

} // namespace
