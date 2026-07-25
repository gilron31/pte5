// Naive schoolbook big-integer arithmetic on the GPU.
//
// One thread handles one operand pair; the batch supplies all the parallelism.
// Limbs are 32-bit so a partial product plus an accumulator limb plus a carry
// always fits in a uint64_t, which keeps the carry handling obvious and avoids
// needing __umul64hi.
#include "check.cuh"
#include "gpu.h"

#include <stdexcept>
#include <string>

namespace pte::gpu {
namespace {

constexpr int kThreadsPerBlock = 128;

int BlocksFor(int items, int threads) {
  return (items + threads - 1) / threads;
}

void ValidateBatchArgs(int nlimbs, int batch) {
  if (nlimbs <= 0 || nlimbs > kMaxLimbs) {
    throw std::invalid_argument("nlimbs must be in [1, " +
                                std::to_string(kMaxLimbs) + "], got " +
                                std::to_string(nlimbs));
  }
  if (batch < 0) {
    throw std::invalid_argument("batch must be non-negative, got " +
                                std::to_string(batch));
  }
}

__global__ void MulBatchKernel(const std::uint32_t *a, const std::uint32_t *b,
                               std::uint32_t *out, int nlimbs, int batch) {
  const int t = blockIdx.x * blockDim.x + threadIdx.x;
  if (t >= batch)
    return;

  const std::uint32_t *x = a + static_cast<long>(t) * nlimbs;
  const std::uint32_t *y = b + static_cast<long>(t) * nlimbs;

  // Per-thread accumulator. Sized by the compile-time bound but only the first
  // 2*nlimbs entries are used.
  std::uint32_t acc[2 * kMaxLimbs];
  for (int i = 0; i < 2 * nlimbs; ++i)
    acc[i] = 0;

  for (int i = 0; i < nlimbs; ++i) {
    std::uint32_t carry = 0;
    for (int j = 0; j < nlimbs; ++j) {
      const std::uint64_t p =
          static_cast<std::uint64_t>(x[i]) * y[j] + acc[i + j] + carry;
      acc[i + j] = static_cast<std::uint32_t>(p);
      carry = static_cast<std::uint32_t>(p >> 32);
    }
    // Outer iteration i is the first to touch position i + nlimbs (iteration i'
    // reaches only up to i' + nlimbs), so the carry lands in a zero limb and
    // needs no further propagation.
    acc[i + nlimbs] = carry;
  }

  std::uint32_t *dst = out + static_cast<long>(t) * 2 * nlimbs;
  for (int i = 0; i < 2 * nlimbs; ++i)
    dst[i] = acc[i];
}

__global__ void AddBatchKernel(const std::uint32_t *a, const std::uint32_t *b,
                               std::uint32_t *out, int nlimbs, int batch) {
  const int t = blockIdx.x * blockDim.x + threadIdx.x;
  if (t >= batch)
    return;

  const std::uint32_t *x = a + static_cast<long>(t) * nlimbs;
  const std::uint32_t *y = b + static_cast<long>(t) * nlimbs;
  std::uint32_t *dst = out + static_cast<long>(t) * (nlimbs + 1);

  std::uint32_t carry = 0;
  for (int i = 0; i < nlimbs; ++i) {
    const std::uint64_t s = static_cast<std::uint64_t>(x[i]) + y[i] + carry;
    dst[i] = static_cast<std::uint32_t>(s);
    carry = static_cast<std::uint32_t>(s >> 32);
  }
  dst[nlimbs] = carry;
}

} // namespace

void mul_batch(const std::uint32_t *a, const std::uint32_t *b,
               std::uint32_t *out, int nlimbs, int batch) {
  ValidateBatchArgs(nlimbs, batch);
  if (batch == 0)
    return;

  const std::size_t in_count = static_cast<std::size_t>(batch) * nlimbs;
  const std::size_t out_count = in_count * 2;

  DeviceBuffer<std::uint32_t> d_a(in_count);
  DeviceBuffer<std::uint32_t> d_b(in_count);
  DeviceBuffer<std::uint32_t> d_out(out_count);
  d_a.CopyFromHost(a);
  d_b.CopyFromHost(b);

  MulBatchKernel<<<BlocksFor(batch, kThreadsPerBlock), kThreadsPerBlock>>>(
      d_a.get(), d_b.get(), d_out.get(), nlimbs, batch);
  PTE_CUDA_CHECK(cudaGetLastError());
  d_out.CopyToHost(out); // synchronizes, and surfaces execution errors
}

void add_batch(const std::uint32_t *a, const std::uint32_t *b,
               std::uint32_t *out, int nlimbs, int batch) {
  ValidateBatchArgs(nlimbs, batch);
  if (batch == 0)
    return;

  const std::size_t in_count = static_cast<std::size_t>(batch) * nlimbs;
  const std::size_t out_count = static_cast<std::size_t>(batch) * (nlimbs + 1);

  DeviceBuffer<std::uint32_t> d_a(in_count);
  DeviceBuffer<std::uint32_t> d_b(in_count);
  DeviceBuffer<std::uint32_t> d_out(out_count);
  d_a.CopyFromHost(a);
  d_b.CopyFromHost(b);

  AddBatchKernel<<<BlocksFor(batch, kThreadsPerBlock), kThreadsPerBlock>>>(
      d_a.get(), d_b.get(), d_out.get(), nlimbs, batch);
  PTE_CUDA_CHECK(cudaGetLastError());
  d_out.CopyToHost(out);
}

int device_count() {
  int count = 0;
  if (cudaGetDeviceCount(&count) != cudaSuccess) {
    cudaGetLastError(); // clear the sticky error so later calls are not
                        // poisoned
    return 0;
  }
  return count;
}

std::string device_name(int index) {
  cudaDeviceProp prop{};
  PTE_CUDA_CHECK(cudaGetDeviceProperties(&prop, index));
  return prop.name;
}

} // namespace pte::gpu
