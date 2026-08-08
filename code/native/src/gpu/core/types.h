// Names for the quantities the buffer types traffic in.
//
// Scope: these are *infrastructure* types -- how big a buffer is, how many
// items a kernel appended. Anything meaning something about the problem (nu
// values, pair keys) belongs to the app that owns it, e.g. `nu_t` in
// gpu/meet/meet.h and `trunc_key_t` inside the meet kernel.
//
// ---------------------------------------------------------------------------
// The file-extension rule for src/gpu, since this header is where a reader
// looks for it:
//
//   .cu / .cuh  iff the file contains *device code* -- __device__/__global__
//               definitions, `<<< >>>`, or device intrinsics like atomicAdd.
//   .cc / .h    everything else, including every call to the CUDA runtime API
//               (cudaMalloc, cudaMemcpy, cudaHostAlloc, ...). That API is a C
//               library; it needs <cuda_runtime.h> and -lcudart, not nvcc.
//
// The rule is worth stating because it is load-bearing: one __device__ function
// in a header promotes *every* translation unit that includes it, transitively,
// to nvcc. That is why DeviceAppendView's data lives in a .h and its push() in
// the matching .cuh.
// ---------------------------------------------------------------------------
//
// A note on `count_t`. It is spelled `unsigned long long` rather than
// std::uint64_t because that is the exact type atomicAdd's 64-bit overload
// takes; on LP64 std::uint64_t is `unsigned long`, a distinct type, and using
// it would either fail to compile or drag in a pointer cast at every atomic.
#pragma once

#include <cstdint>

namespace pte::gpu {

// Indices and sizes of batched work: CSR offsets, norm counts, buffer
// capacities. 64-bit because a full campaign batch has far more than 2^31
// points, and a 32-bit overflow here would silently truncate a result rather
// than fail.
using index_t = std::int64_t;

// Counter incremented by atomicAdd. See the header comment for the spelling.
using count_t = unsigned long long;

} // namespace pte::gpu
