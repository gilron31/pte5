// Host-side API for the CUDA kernels.
//
// This header is deliberately plain C++: no <cuda_runtime.h>, no __global__, no
// CUDA types in any signature. That keeps it includable from ordinary g++ TUs
// (the demo app, the Python bindings) and keeps the GPU unit tests free of any
// dependency on the CPU library.
//
// Big-integer convention, shared with src/cpu/bigint.h: an operand is a
// fixed-width array of 32-bit limbs, least-significant limb first.
#ifndef PTE_GPU_GPU_H_
#define PTE_GPU_GPU_H_

#include <cstdint>
#include <string>

namespace pte::gpu {

// Largest operand width the kernels accept, in 32-bit limbs (1024 bits). The
// schoolbook kernel keeps its accumulator in a per-thread local array, so the
// bound has to be a compile-time constant.
inline constexpr int kMaxLimbs = 32;

// Number of usable CUDA devices; 0 if the runtime reports none or errors out.
// Never throws, so tests can use it to decide whether to skip.
int device_count();

// Name of device `index`. Throws std::runtime_error if unavailable.
std::string device_name(int index = 0);

// Batched big-integer multiply: for each t < batch,
//   out[t*2*nlimbs ...] = a[t*nlimbs ...] * b[t*nlimbs ...]
// Naive schoolbook, one thread per pair. `out` needs batch*2*nlimbs limbs.
// Throws std::invalid_argument unless 0 < nlimbs <= kMaxLimbs and batch >= 0.
void mul_batch(const std::uint32_t *a, const std::uint32_t *b,
               std::uint32_t *out, int nlimbs, int batch);

// Batched big-integer add with carry-out: for each t < batch,
//   out[t*(nlimbs+1) ...] = a[t*nlimbs ...] + b[t*nlimbs ...]
// `out` needs batch*(nlimbs+1) limbs; the extra limb holds the final carry.
void add_batch(const std::uint32_t *a, const std::uint32_t *b,
               std::uint32_t *out, int nlimbs, int batch);

// Row-major (M x K) * (K x N) = (M x N), double precision, tiled kernel.
// C needs M*N doubles.
void matmul(const double *A, const double *B, double *C, int M, int N, int K);

} // namespace pte::gpu

#endif // PTE_GPU_GPU_H_
