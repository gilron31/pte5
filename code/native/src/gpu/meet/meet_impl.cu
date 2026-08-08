// The E_4 meet kernel.
//
// Scope, which the name is meant to carry:
//
//   "E4"     - it looks for *quadruples*. An E_4 solution is four points, so a
//              4-subset of a norm's point set is exactly one candidate. This
//              kernel says nothing about E_5, which needs eight points and a
//              different search entirely.
//   "Direct" - it compares every 4-subset directly rather than hashing pair
//              keys. That is the right choice only while the point-set size D
//              is small: the work is C(D,4)*3 ~ D^4/8, which is 3 comparisons
//              at D = 4 (about three quarters of admissible norms) but 1.6e8
//              at the D = 192 that B = 1e12 reaches. A hashing variant is
//              needed for large D; see drafts/gpu_algorithms.md A4.
//
// One thread per norm. Keys are compared modulo 2^64 by unsigned wraparound;
// equal keys are always congruent, so there are no false negatives, only false
// positives for the host to reject in wider arithmetic.
//
// No disjointness test appears anywhere: nu is injective on a point set (see
// src/cpu/gaussian/gaussian.h), so two distinct pairs with equal key cannot
// share a point.
#include "gpu/core/device_append_buffer.h"
#include "gpu/core/launch.cuh"
#include "gpu/meet/meet_impl.h"

namespace pte::gpu {
namespace {

// The pair key nu_a^2 + nu_b^2 reduced modulo 2^64 by unsigned wraparound.
//
// The name says "truncated" because that is the only thing standing between
// this kernel and a wrong answer: above pte::cpu::kExactKeyBound the exact key
// needs more than 64 bits, so a match here is a *candidate*, and the host must
// re-check it in wider arithmetic. Never let this type escape the kernel under
// a name that suggests it is the key.
using trunc_key_t = std::uint64_t;

__device__ inline trunc_key_t PairKey(nu_t x, nu_t y) {
  return static_cast<trunc_key_t>(x) * x + static_cast<trunc_key_t>(y) * y;
}

template <int kThreads, int kMinBlocksPerSM>
__global__ __launch_bounds__(kThreads, kMinBlocksPerSM) void E4MeetDirectKernel(
    const nu_t *__restrict__ nu, const index_t *__restrict__ offsets,
    index_t num_norms, DeviceAppendView<MeetHit> out) {
  const index_t b = static_cast<index_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (b >= num_norms)
    return;

  const index_t begin = offsets[b];
  const int n = static_cast<int>(offsets[b + 1] - begin);
  const nu_t *__restrict__ v = nu + begin;

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const trunc_key_t kij = PairKey(v[i], v[j]);
      for (int k = j + 1; k < n; ++k) {
        const trunc_key_t kik = PairKey(v[i], v[k]);
        for (int l = k + 1; l < n; ++l) {
          // The three pairings of {i,j,k,l}, tested independently: a 4-subset
          // can match in more than one pairing, and the host reference reports
          // each, so the two must agree.
          if (kij == PairKey(v[k], v[l]))
            out.push(MeetHit{b, i, j, k, l});
          if (kik == PairKey(v[j], v[l]))
            out.push(MeetHit{b, i, k, j, l});
          if (PairKey(v[i], v[l]) == PairKey(v[j], v[k]))
            out.push(MeetHit{b, i, l, j, k});
        }
      }
    }
  }
}

} // namespace

void LaunchE4MeetDirect(const DeviceBuffer<nu_t> &d_nu,
                        const DeviceBuffer<index_t> &d_offsets,
                        DeviceAppendView<MeetHit> out) {
  if (d_offsets.count() < 2)
    return; // fewer than two offsets describes no norms
  const auto num_norms = static_cast<index_t>(d_offsets.count() - 1);
  const index_t blocks =
      (num_norms + kMeetThreadsPerBlock - 1) / kMeetThreadsPerBlock;
  // The one place raw device pointers appear: a kernel argument list cannot
  // take anything else.
  PTE_CUDA_LAUNCH_KERNEL(
      (E4MeetDirectKernel<kMeetThreadsPerBlock, kMeetMinBlocksPerSM>),
      static_cast<unsigned>(blocks), kMeetThreadsPerBlock, d_nu.get(),
      d_offsets.get(), num_norms, out);
}

} // namespace pte::gpu
