// Device-side contract for the champion meet.
//
// Plain C++ despite naming device concepts: this declares a host function that
// launches a kernel, and DeviceAppendView is a POD. Only champion_impl.cu needs
// nvcc.
#pragma once

#include "core/uint.h"
#include "gpu/core/device_append_buffer.h"
#include "gpu/core/device_buffer.h"
#include "gpu/core/types.h"

#include <cstdint>

namespace pte::gpu {

// A pair of pairs the device believes have equal keys, as indices into the
// *sorted* point set. A candidate, not a result: the signature the device
// compares is narrower than the key (see champion.h), so the host re-checks
// each one in exact 256-bit arithmetic before calling it a hit.
struct ChampionCandidate {
  std::uint32_t a = 0, b = 0, c = 0, d = 0;
};

// Points per slot-index field, and therefore the ceiling on D. 2^20 against the
// 2^19 a 20-prime champion produces, which is the largest point set anyone has
// proposed building.
inline constexpr int kChampionIndexBits = 20;
inline constexpr std::uint32_t kChampionMaxPoints = 1u << kChampionIndexBits;

// One warp per point i, so a warp's 32 lanes read 32 consecutive sq[j] -- 1 KB
// contiguous, fully coalesced -- and every lane in a warp shares one i and one
// j-range, which is where the divergence would otherwise be.
//
// A block per i was the obvious alternative and is worse: the j-runs average
// ~64 pairs at the sized window, so a 256-thread block would leave three
// quarters of its threads idle.
inline constexpr int kChampionThreadsPerBlock = 256;

// Runs one value window over buffers that already live on the device.
//
// `d_table` must be zeroed and `out` reset by the caller before each call --
// the kernel appends, it does not clear. `d_overflow` counts probe sequences
// that ran to the limit; a nonzero value means the table saturated and the
// window's result is incomplete, which the caller must treat as an error rather
// than as "no collisions here".
//
// Takes the owning buffers rather than raw pointers, so the point count is
// d_sq.count() and the slot count is d_table.count() instead of parameters a
// caller could contradict (REVIEW.md D17). Raw device pointers appear exactly
// once, in the kernel's argument list.
void LaunchChampionWindow(const DeviceBuffer<core::U256> &d_sq,
                          const DeviceBuffer<std::uint32_t> &d_lo,
                          const DeviceBuffer<std::uint32_t> &d_hi,
                          DeviceBuffer<count_t> &d_table,
                          DeviceBuffer<count_t> &d_overflow,
                          DeviceAppendView<ChampionCandidate> out);

} // namespace pte::gpu
