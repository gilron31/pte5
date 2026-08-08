// CUDA error checking.
//
// Only included from .cu translation units. The owning buffer types live in
// device_buffer.cuh and host_buffer.cuh.
#pragma once

#include <cuda_runtime.h>

#include <stdexcept>
#include <string>
// GILRO: this file has no reason to be a .cuh, it should be an .h file.
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

// Launch a kernel and check the launch immediately.
//
// A launch is asynchronous, so cudaGetLastError() right after it reports
// configuration errors (bad grid or block shape, too much shared memory) but
// not execution errors; those surface at the next synchronizing call, which is
// why every copy back is PTE_CUDA_CHECK'd too. Going through this macro means
// no launch can silently skip the configuration check.
#define PTE_CUDA_LAUNCH_KERNEL(kernel, blocks, threads, ...)                   \
  do {                                                                         \
    (kernel)<<<(blocks), (threads)>>>(__VA_ARGS__);                            \
    PTE_CUDA_CHECK(cudaGetLastError());                                        \
  } while (0)
