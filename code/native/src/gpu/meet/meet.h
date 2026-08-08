// The E_4 meet: search S-A of drafts/gpu_algorithms.md.
//
// Plain C++ host API: no <cuda_runtime.h>, no CUDA types. The implementation
// (meet.cu) owns the staging and device buffers and calls into meet_impl.cuh.
#pragma once

#include "gpu/core/types.h"

#include <cstdint>

namespace pte::gpu {

// One found configuration: within norm `norm`, the point indices split as
// {a,b} | {c,d}. Indices are local to that norm's slice of the nu array.
//
// A quadruple is a candidate, not a confirmed solution: keys are compared
// modulo 2^64, so the host must re-check in wider arithmetic (see below).
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
// Layout: norm i owns nu[offsets[i] .. offsets[i+1]). `offsets` therefore has
// num_norms + 1 entries, must be non-decreasing, and must end at nu_count.
//
// Arithmetic: keys are compared modulo 2^64 by unsigned wraparound. That is
// exact while the caller's nu values stay inside pte::cpu::kExactKeyBound and a
// sound filter beyond it -- equal keys are always congruent, so there are no
// false negatives, only false positives that the host must reject in wider
// arithmetic.
//
// Returns the total number of hits found, which may exceed max_hits; in that
// case the first max_hits were written and the caller should retry with a
// buffer of the returned size. Throws std::invalid_argument on a negative
// size, a null pointer where one is required, or a malformed offsets array.
index_t meet_batch(const nu_t *nu, index_t nu_count, const index_t *offsets,
                   index_t num_norms, MeetHit *hits, index_t max_hits);

} // namespace pte::gpu
