// CUDA error checking.
//
// Ordinary C++, compiled by the host compiler: the CUDA *runtime API* is a C
// library, so calling it needs <cuda_runtime.h> and -lcudart, not nvcc. Only
// device code needs nvcc. See the extension rule at the top of
// gpu/core/types.h.
//
// The kernel-launch macro lives in launch.cuh, because `<<< >>>` is device
// syntax and would pull every includer into nvcc.
#pragma once

#include <cuda_runtime.h>

#include <stdexcept>
#include <string>

namespace pte::gpu {

inline void ThrowIfCudaError(cudaError_t err, const char *expr,
                             const char *file, int line) {
  if (err != cudaSuccess) {
    throw std::runtime_error(std::string(file) + ":" + std::to_string(line) +
                             ": " + expr + " failed: " + cudaGetErrorName(err) +
                             " (" + cudaGetErrorString(err) + ")");
  }
}

} // namespace pte::gpu

// Wrap every runtime call.
#define PTE_CUDA_CHECK(expr)                                                   \
  ::pte::gpu::ThrowIfCudaError((expr), #expr, __FILE__, __LINE__)
