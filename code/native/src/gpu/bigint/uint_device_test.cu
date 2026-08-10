// core/uint.h computed on the device, checked against core/uint.h computed on
// the host.
//
// This is the obligation core/hostdev.h attaches to sharing a type across the
// seam. `add` and `sub` have *two bodies* -- an inline PTX carry chain under
// __CUDA_ARCH__ and a portable uint64 accumulate otherwise -- and the whole
// arrangement is sound only while the two agree on every input. Nothing in the
// language enforces that; this file does.
//
// It lives in gpu/bigint because that is the app the plan promotes for S-E
// ("gpu/bigint/ -- Promote: S-E needs a fixed 4-limb type built on this
// convention"), and because the check is device-side by nature: it must launch
// a kernel, so it is a .cu, and the host oracle it compares against is the same
// header compiled for the host -- which is exactly the comparison wanted.
//
// The values are not random-uniform. A carry chain is wrong at carries, so the
// corpus is built to produce them: all-ones limbs, single-bit limbs, and pairs
// chosen to carry out of limb k for every k.
#include "core/uint.h"

#include "gpu/core/check.h"
#include "gpu/core/device.h"
#include "gpu/core/device_buffer.h"
#include "gpu/core/launch.cuh"

#include <gtest/gtest.h>

#include <cstdint>
#include <random>
#include <vector>

namespace pte::core {
namespace {

// One operand pair per thread; every operation applied to it.
struct Result {
  U256 sum;
  U256 diff;
  U256 product; // low 256 bits of a 128x128 product of the low halves
  int cmp;
  int bits;
};

__global__ void ApplyAll(const U256 *__restrict__ a, const U256 *__restrict__ b,
                         Result *__restrict__ out, int n) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n)
    return;
  const U256 x = a[i], y = b[i];
  Result r;
  r.sum = add(x, y);
  r.diff = sub(x, y);
  U128 xl, yl;
  for (int k = 0; k < 4; ++k) {
    xl.w[k] = x.w[k];
    yl.w[k] = y.w[k];
  }
  r.product = mul(xl, yl);
  r.cmp = cmp(x, y);
  r.bits = bit_length(x);
  out[i] = r;
}

std::vector<U256> Operands() {
  std::vector<U256> v;
  const std::uint32_t m = 0xffffffffu;
  // Carry-producing shapes: k limbs of all ones, then a 1 that must ripple
  // through every one of them.
  for (int k = 0; k <= 8; ++k) {
    U256 x{};
    for (int i = 0; i < k; ++i)
      x.w[i] = m;
    v.push_back(x);
    U256 y{};
    if (k < 8)
      y.w[k] = 1;
    else
      y.w[0] = 1;
    v.push_back(y);
  }
  std::mt19937_64 rng(90210);
  for (int k = 0; k < 200; ++k) {
    U256 x{};
    const int limbs = 1 + (k % 8);
    for (int i = 0; i < limbs; ++i)
      x.w[i] = static_cast<std::uint32_t>(rng());
    v.push_back(x);
  }
  return v;
}

TEST(UIntDevice, DeviceAndHostAgreeOnEveryOperation) {
  if (pte::gpu::device_count() == 0)
    GTEST_SKIP() << "no CUDA device";

  const std::vector<U256> corpus = Operands();
  std::vector<U256> a, b;
  for (std::size_t i = 0; i < corpus.size(); ++i)
    for (std::size_t j = 0; j < corpus.size(); j += 7) {
      a.push_back(corpus[i]);
      b.push_back(corpus[j]);
    }
  const int n = static_cast<int>(a.size());
  ASSERT_GT(n, 1000);

  pte::gpu::DeviceBuffer<U256> d_a(a.size()), d_b(b.size());
  pte::gpu::DeviceBuffer<Result> d_out(a.size());
  d_a.CopyFromHost(a.data());
  d_b.CopyFromHost(b.data());

  constexpr int kThreads = 128;
  PTE_CUDA_LAUNCH_KERNEL(ApplyAll, (n + kThreads - 1) / kThreads, kThreads,
                         d_a.get(), d_b.get(), d_out.get(), n);
  std::vector<Result> got(a.size());
  d_out.CopyToHost(got.data());

  int mismatches = 0;
  for (int i = 0; i < n && mismatches < 5; ++i) {
    const U256 &x = a[static_cast<std::size_t>(i)];
    const U256 &y = b[static_cast<std::size_t>(i)];
    const Result &r = got[static_cast<std::size_t>(i)];
    U128 xl, yl;
    for (int k = 0; k < 4; ++k) {
      xl.w[k] = x.w[k];
      yl.w[k] = y.w[k];
    }
    if (!(r.sum == add(x, y))) {
      ADD_FAILURE() << "device add disagrees at i=" << i;
      ++mismatches;
    }
    if (!(r.diff == sub(x, y))) {
      ADD_FAILURE() << "device sub disagrees at i=" << i;
      ++mismatches;
    }
    if (!(r.product == mul(xl, yl))) {
      ADD_FAILURE() << "device mul disagrees at i=" << i;
      ++mismatches;
    }
    if (r.cmp != cmp(x, y)) {
      ADD_FAILURE() << "device cmp disagrees at i=" << i;
      ++mismatches;
    }
    if (r.bits != bit_length(x)) {
      ADD_FAILURE() << "device bit_length disagrees at i=" << i;
      ++mismatches;
    }
  }
}

// The layout half of the contract: a U256 written by the host and read back
// through device code must be the same 32 bytes. A padding or alignment
// surprise here would corrupt every key silently rather than fail to compile.
__global__ void RoundTrip(const U256 *in, U256 *out, int n) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    out[i] = in[i];
}

TEST(UIntDevice, LayoutSurvivesTheRoundTrip) {
  if (pte::gpu::device_count() == 0)
    GTEST_SKIP() << "no CUDA device";
  static_assert(sizeof(U256) == 32, "");

  const std::vector<U256> in = Operands();
  const int n = static_cast<int>(in.size());
  pte::gpu::DeviceBuffer<U256> d_in(in.size()), d_out(in.size());
  d_in.CopyFromHost(in.data());
  constexpr int kThreads = 128;
  PTE_CUDA_LAUNCH_KERNEL(RoundTrip, (n + kThreads - 1) / kThreads, kThreads,
                         d_in.get(), d_out.get(), n);
  std::vector<U256> out(in.size());
  d_out.CopyToHost(out.data());
  for (std::size_t i = 0; i < in.size(); ++i)
    EXPECT_TRUE(in[i] == out[i]) << "limb layout changed crossing the bus";
}

} // namespace
} // namespace pte::core
