// The E_4 "meet": for each norm, find four points whose nu values satisfy
//     nu[i]^2 + nu[j]^2 == nu[k]^2 + nu[l]^2
// over disjoint pairs {i,j} and {k,l}.
//
// One thread per norm. Each walks the 4-subsets of its point set and tests the
// three ways to split each into two pairs. This is the direct-comparison form,
// which is the right one here: the point-set size D is 4 for about three
// quarters of admissible norms and at most a few hundred at the bounds we
// enumerate, so a hash table would cost more than it saves.
//
// No disjointness test appears anywhere. nu is injective on a point set (see
// src/cpu/gaussian.h), so two distinct pairs with equal key cannot share a
// point.
//
// Keys are compared modulo 2^64 by unsigned wraparound. Equal keys are always
// congruent, so the kernel has no false negatives -- only false positives, left
// for the host to reject in wider arithmetic.
#include "check.cuh"
#include "gpu.h"

#include <stdexcept>

namespace pte::gpu {
namespace {

__device__ inline unsigned long long PairKey(unsigned long long x,
                                             unsigned long long y) {
  return x * x + y * y; // wraps; see the file comment
}

// Appends (norm, a, b, c, d) when `matched`. A slot is reserved even once the
// buffer is full, so the returned count tells the caller how much room a retry
// needs.
__device__ inline void Emit(bool matched, int norm, int a, int b, int c, int d,
                            int *hits, int max_hits, unsigned int *hit_count) {
  if (!matched)
    return;
  const unsigned int slot = atomicAdd(hit_count, 1u);
  if (slot < static_cast<unsigned int>(max_hits)) {
    int *out = hits + static_cast<std::size_t>(slot) * kMeetHitStride;
    out[0] = norm;
    out[1] = a;
    out[2] = b;
    out[3] = c;
    out[4] = d;
  }
}

__global__ void MeetKernel(const unsigned long long *nu, const int *offsets,
                           int num_norms, int *hits, int max_hits,
                           unsigned int *hit_count) {
  const int b = blockIdx.x * blockDim.x + threadIdx.x;
  if (b >= num_norms)
    return;

  const int begin = offsets[b];
  const int n = offsets[b + 1] - begin;
  const unsigned long long *v = nu + begin;

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const unsigned long long kij = PairKey(v[i], v[j]);
      for (int k = j + 1; k < n; ++k) {
        const unsigned long long kik = PairKey(v[i], v[k]);
        for (int l = k + 1; l < n; ++l) {
          // The three pairings of {i,j,k,l}, tested independently: a 4-subset
          // can match in more than one pairing, and the host reference reports
          // each, so the two must agree.
          Emit(kij == PairKey(v[k], v[l]), b, i, j, k, l, hits, max_hits,
               hit_count);
          Emit(kik == PairKey(v[j], v[l]), b, i, k, j, l, hits, max_hits,
               hit_count);
          Emit(PairKey(v[i], v[l]) == PairKey(v[j], v[k]), b, i, l, j, k, hits,
               max_hits, hit_count);
        }
      }
    }
  }
}

} // namespace

int meet_batch(const std::uint64_t *nu, const std::int32_t *offsets,
               int num_norms, std::int32_t *hits, int max_hits) {
  if (num_norms < 0 || max_hits < 0)
    throw std::invalid_argument("meet_batch: negative size");
  if (num_norms == 0)
    return 0;
  if (offsets == nullptr || nu == nullptr)
    throw std::invalid_argument("meet_batch: null input");
  if (max_hits > 0 && hits == nullptr)
    throw std::invalid_argument("meet_batch: null hits buffer");

  const std::size_t nu_count = static_cast<std::size_t>(offsets[num_norms]);

  DeviceBuffer<unsigned long long> d_nu(nu_count);
  DeviceBuffer<int> d_offsets(static_cast<std::size_t>(num_norms) + 1);
  DeviceBuffer<int> d_hits(static_cast<std::size_t>(max_hits) * kMeetHitStride);
  DeviceBuffer<unsigned int> d_count(1);

  if (nu_count > 0)
    d_nu.CopyFromHost(reinterpret_cast<const unsigned long long *>(nu));
  d_offsets.CopyFromHost(reinterpret_cast<const int *>(offsets));
  d_count.Zero();
  d_hits.Zero(); // unwritten tail is copied back too; keep it deterministic

  constexpr int kThreads = 128;
  const int blocks = (num_norms + kThreads - 1) / kThreads;
  MeetKernel<<<blocks, kThreads>>>(d_nu.get(), d_offsets.get(), num_norms,
                                   d_hits.get(), max_hits, d_count.get());
  PTE_CUDA_CHECK(cudaGetLastError());

  unsigned int found = 0;
  d_count.CopyToHost(&found); // synchronizes, and surfaces execution errors
  const int written =
      static_cast<int>(found) < max_hits ? static_cast<int>(found) : max_hits;
  if (written > 0)
    d_hits.CopyToHost(reinterpret_cast<int *>(hits));
  return static_cast<int>(found);
}

} // namespace pte::gpu
