// Device discovery. Plain C++: no <cuda_runtime.h>, no CUDA types.
#pragma once

#include <string>

namespace pte::gpu {

// Number of usable CUDA devices; 0 if the runtime reports none or errors out.
// Never throws, so tests can use it to decide whether to skip.
int device_count();

// Name of device `index`. Throws std::runtime_error if unavailable.
std::string device_name(int index = 0);

} // namespace pte::gpu
