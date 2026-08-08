// Kernel launch.
//
// Separate from check.h because the `<<< >>>` in the macro body is device
// syntax: any translation unit that *expands* this macro must be compiled by
// nvcc. Keeping it here means including check.h costs a host-only file nothing.
#pragma once

#include "gpu/core/check.h"

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
