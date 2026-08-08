// The E_4 meet: search S-A of drafts/gpu_algorithms.md.
//
// Plain C++ host API: no CUDA types in any signature. The implementation
// (meet.cc) owns the staging and device buffers and calls into meet_impl.h.
#pragma once

#include "gpu/core/types.h"

#include <cstdint>
#include <span>

namespace pte::gpu {

// The point invariant nu(z) = Re(z^2) = m^2 - n^2, non-negative on canonized
// points and bounded by L1. Owned by the meet rather than by gpu/core: it means
// something about *this* problem, and nothing else in src/gpu traffics in it.
//
// (The truncated pair key nu_a^2 + nu_b^2 mod 2^64 is narrower still -- it
// never crosses this API and is declared inside meet_impl.cu.)
using nu_t = std::uint64_t;

// One found configuration: within norm `norm`, the point indices split as
// {a,b} | {c,d}. Indices are local to that norm's slice of the nu array.
//
// A quadruple is a candidate, not a confirmed solution: keys are compared
// modulo 2^64, so the host must re-check in wider arithmetic (see below).
//
// Deliberately not pte::cpu::Quadruple. Reusing that type would make
// pte_gpu_meet link pte_cpu, and the CPU enumerator is the oracle this kernel
// is checked against -- the two halves must be able to build and be tested
// independently. See REVIEW.md D3.
struct MeetHit {
  index_t norm = 0;
  std::int32_t a = 0;
  std::int32_t b = 0;
  std::int32_t c = 0;
  std::int32_t d = 0;

  friend bool operator==(const MeetHit &x, const MeetHit &y) {
    return x.norm == y.norm && x.a == y.a && x.b == y.b && x.c == y.c &&
           x.d == y.d;
  }
};

// Finds, for every norm in the batch, all index quadruples with
//
//     nu[a]^2 + nu[b]^2 == nu[c]^2 + nu[d]^2
//
// over disjoint pairs {a,b} and {c,d}. No disjointness test is performed or
// needed: nu is injective on a point set, so two distinct pairs with equal key
// cannot share a point.
//
// Layout: norm i owns nu[offsets[i] .. offsets[i+1]), so `offsets` has one more
// entry than there are norms, must be non-decreasing, and must end at or before
// nu.size(). A batch of n norms therefore passes n+1 offsets; an empty or
// single-element `offsets` describes no work and returns 0.
//
// Arithmetic: keys are compared modulo 2^64 by unsigned wraparound. That is
// exact while the caller's nu values stay inside pte::cpu::kExactKeyBound and a
// sound filter beyond it -- equal keys are always congruent, so there are no
// false negatives, only false positives that the host must reject in wider
// arithmetic.
//
// Returns the total number of hits found, which may exceed hits.size(); in that
// case the first hits.size() were written and the caller should retry with a
// buffer of the returned size. Throws std::invalid_argument on a malformed
// offsets array.
//
// Spans rather than (pointer, count) triples: the count of a span cannot
// disagree with its data, which is the entire content of the validation that
// used to live at the top of the implementation. (The device-side contract in
// meet_impl.h keeps raw pointers -- those address device memory, and std::span
// has no __device__ members.)
index_t meet_batch(std::span<const nu_t> nu, std::span<const index_t> offsets,
                   std::span<MeetHit> hits);

} // namespace pte::gpu
