// Device-side contract for the E_4 meet.
//
// Plain C++ despite naming a device concept: this declares a host function that
// happens to launch a kernel, and DeviceAppendView is a POD. Only meet_impl.cu
// needs nvcc.
#pragma once

#include "gpu/core/device_append_buffer.h"
#include "gpu/core/device_buffer.h"
#include "gpu/core/types.h"
#include "gpu/meet/meet.h"

namespace pte::gpu {

// Launch configuration for the meet, in one place so the __launch_bounds__ on
// the kernel and the launch itself cannot drift apart.
//
// Both values are template parameters of the kernel, so a bucketed variant can
// be tuned per D without duplicating the body.
inline constexpr int kMeetThreadsPerBlock = 128;

// The second __launch_bounds__ argument, minBlocksPerMultiprocessor, is a
// register-pressure knob, not a scheduling hint: it tells nvcc how many blocks
// must fit concurrently, and nvcc will spill to local memory to keep that
// promise. Measured on sm_89 (65536 registers/SM, 1536 threads/SM) at 128
// threads per block, with -Xptxas -v:
//
//   arg   registers   spilled   blocks/SM   occupancy
//   ----------------------------------------------------
//   none     56          0          9          75%
//    1       66          0          7          58%
//    6       62          0          8          67%
//    9       56          0          9          75%
//   12       40         68 B       12         100%
//
// Note that 1 is *not* "no constraint": it promises only one resident block, so
// nvcc feels free to spend 66 registers and occupancy drops. 12 buys full
// occupancy by spilling, which for a loop this ALU-heavy is the wrong trade.
//
// 9 reproduces what nvcc picks unaided, and pinning it means a later edit to
// the kernel body cannot silently drift to the 7-block configuration. It is
// tuned for sm_89, the one architecture this project targets; re-measure before
// changing CMAKE_CUDA_ARCHITECTURES, and re-measure per bucket when the
// D-bucketed kernels land.
//
// Occupancy is a proxy here, not a timing: the kernel is microseconds against a
// transfer-bound demo, so wall-clock differences between these rows are not
// measurable yet.
inline constexpr int kMeetMinBlocksPerSM = 9;

// Launches the meet over buffers that already live on the device, appending
// every candidate to `out`. Does not synchronize; the caller's copy back does.
//
// Takes the owning buffers, not raw device pointers: the buffer carries its own
// extent, so the norm count is `d_offsets.count() - 1` rather than a fourth
// parameter a caller could contradict -- the same argument as the std::span
// signature on meet_batch. Raw pointers appear exactly once, inside the launch.
void LaunchE4MeetDirect(const DeviceBuffer<nu_t> &d_nu,
                        const DeviceBuffer<index_t> &d_offsets,
                        DeviceAppendView<MeetHit> out);

} // namespace pte::gpu
