#include "gpu/core/device.h"

#include "gpu/core/check.h"

namespace pte::gpu {

int device_count() {
  int count = 0;
  // Deliberately not PTE_CUDA_CHECK'd: a machine with no driver at all reports
  // an error here, and callers use this to decide whether to skip.
  if (cudaGetDeviceCount(&count) != cudaSuccess) {
    cudaGetLastError(); // clear the sticky error so later calls are not
                        // poisoned
    return 0;
  }
  return count;
}

std::string device_name(int index) {
  cudaDeviceProp prop{};
  PTE_CUDA_CHECK(cudaGetDeviceProperties(&prop, index));
  return prop.name;
}

} // namespace pte::gpu
