// Device-side contract for the E_4 meet. Only included from .cu files.
#pragma once

#include "gpu/core/device_append_buffer.cuh"
#include "gpu/core/types.h"
#include "gpu/meet/meet.h"

namespace pte::gpu {

// Threads per block for the meet. Exposed so the __launch_bounds__ on the
// kernel and the launch configuration cannot drift apart.
inline constexpr int kMeetThreadsPerBlock = 128;

// Launches the meet over buffers that already live on the device, appending
// every candidate to `out`. Does not synchronize; the caller's copy back does.
void LaunchE4MeetDirect(const nu_t *d_nu, const index_t *d_offsets,
                        index_t num_norms, DeviceAppendView<MeetHit> out);

} // namespace pte::gpu

// GILRO: SHould be a .h