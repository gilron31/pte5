// Names for the quantities that cross the host/device boundary.
//
// Plain C++: includable from ordinary g++ translation units as well as from
// .cu files, so host and device agree on one spelling per quantity and no
// reinterpret_cast is needed anywhere.
//
// A note on `count_t`. It is spelled `unsigned long long` rather than
// std::uint64_t because that is the exact type atomicAdd's 64-bit overload
// takes; on LP64 std::uint64_t is `unsigned long`, a distinct type, and using
// it would either fail to compile or drag in a pointer cast at every atomic.
#pragma once

#include <cstdint>

namespace pte::gpu {

// Point invariant nu(z) = Re(z^2), and the pair key nu_a^2 + nu_b^2.
//
// The key is deliberately computed modulo 2^64 by unsigned wraparound: equal
// keys are always congruent, so the comparison has no false negatives, only
// false positives for the host to reject in wider arithmetic.
using nu_t = std::uint64_t;
using key_t = std::uint64_t;

// Indices and sizes of batched work: CSR offsets, norm counts, buffer
// capacities. 64-bit because a full campaign batch has far more than 2^31
// points, and a 32-bit overflow here would silently truncate a result rather
// than fail.
using index_t = std::int64_t;

// Counter incremented by atomicAdd. See the header comment for the spelling.
using count_t = unsigned long long;

} // namespace pte::gpu
