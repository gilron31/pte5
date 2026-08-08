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
// src/cpu/gaussian.h), so two distinct pairs with equal key cannot share a
// point.
#include "gpu/core/check.cuh"
#include "gpu/meet/meet_impl.cuh"

namespace pte::gpu {
namespace {

__device__ inline key_t PairKey(nu_t x, nu_t y) {
  return static_cast<key_t>(x) * x + static_cast<key_t>(y) * y; // wraps
}

// GILRO: lauch bounds also blocks per sm, both as template argument.
__global__ __launch_bounds__(kMeetThreadsPerBlock) void E4MeetDirectKernel(
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
      const key_t kij = PairKey(v[i], v[j]);
      for (int k = j + 1; k < n; ++k) {
        const key_t kik = PairKey(v[i], v[k]);
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

void LaunchE4MeetDirect(const nu_t *d_nu, const index_t *d_offsets,
                        index_t num_norms, DeviceAppendView<MeetHit> out) {
  if (num_norms <= 0)
    return;
  const index_t blocks =
      (num_norms + kMeetThreadsPerBlock - 1) / kMeetThreadsPerBlock;
  PTE_CUDA_LAUNCH_KERNEL(E4MeetDirectKernel, static_cast<unsigned>(blocks),
                         kMeetThreadsPerBlock, d_nu, d_offsets, num_norms, out);
}

} // namespace pte::gpu
