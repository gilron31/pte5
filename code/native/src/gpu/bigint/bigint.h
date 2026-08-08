// Batched big-integer arithmetic on the GPU. Plain C++: no <cuda_runtime.h>,
// no CUDA types.
//
// Operand convention, shared with src/cpu/bigint.h: a big integer is a
// fixed-width array of 32-bit limbs, least-significant limb first.
#pragma once

#include <cstdint>

namespace pte::gpu {

// Largest operand width the kernels accept, in 32-bit limbs (1024 bits). The
// schoolbook kernel keeps its accumulator in a per-thread local array, so the
// bound has to be a compile-time constant.
inline constexpr int kMaxLimbs = 32;

// For each t < batch,
//   out[t*2*nlimbs ...] = a[t*nlimbs ...] * b[t*nlimbs ...]
// Naive schoolbook, one thread per pair. `out` needs batch*2*nlimbs limbs.
// Throws std::invalid_argument unless 0 < nlimbs <= kMaxLimbs and batch >= 0.
void mul_batch(const std::uint32_t *a, const std::uint32_t *b,
               std::uint32_t *out, int nlimbs, int batch);

// For each t < batch,
//   out[t*(nlimbs+1) ...] = a[t*nlimbs ...] + b[t*nlimbs ...]
// `out` needs batch*(nlimbs+1) limbs; the extra limb holds the final carry.
void add_batch(const std::uint32_t *a, const std::uint32_t *b,
               std::uint32_t *out, int nlimbs, int batch);

} // namespace pte::gpu
