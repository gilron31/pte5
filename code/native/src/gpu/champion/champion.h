// S-E, the champion meet, on the device.
//
// Plain C++ host API: no CUDA types in any signature, so the demo, the tests
// and the bindings are all compiled by g++.
//
// ---------------------------------------------------------------------------
// What runs where, and why this shape
// ---------------------------------------------------------------------------
//
// The search is C(D,2) = 1.374e11 pair keys for one champion norm, sliced into
// value windows (core/window.h) that each fit in a device-resident hash table.
// Per window: enumerate the pairs, insert their keys, report equal ones.
//
// The device does the part that is a *bandwidth* problem and nothing else:
//
//   host    sort nu, square it, plan the windows (O(D) per window), confirm
//           every candidate in exact 256-bit arithmetic, map back to the
//           caller's indexing
//   device  for every pair in the window: one 256-bit add, one hash, one
//           atomicCAS probe
//
// bench/gpu P3 and P4 are why it is a hash and not a sort:
//
//   CUB radix sort, 8+8 B pairs   8.4e8 key/s   -> 163 s, and ~2 GB of temp
//   atomicCAS open addressing     1.0e9 ins/s   -> 138 s, table only
//
// The hash wins on both time and memory, and the memory is what actually
// decides it: without CUB's temp storage the window count is free to be chosen
// for convenience rather than to make a sort fit.
//
// ---------------------------------------------------------------------------
// Soundness: filter on the device, arbiter on the host
// ---------------------------------------------------------------------------
//
// This is cpu/core/precision.h's three-layer policy, at champion width. The
// device compares a *signature* of the key -- the hash bits that pick the slot,
// plus a fingerprint stored beside the pair indices -- which is a sound filter
// in exactly one direction:
//
//   equal keys  =>  identical 256-bit words  =>  identical hash
//                =>  same slot sequence and same fingerprint  =>  found.
//
// So there are **no false negatives**. Unequal keys can share a signature, and
// those false positives are re-checked here on the host with full 256-bit
// arithmetic before anything is returned. `ChampionStats::false_positives`
// counts them, so the rate is observable rather than assumed.
//
// Two failure modes are *not* silent, by construction, because both would cost
// a solution:
//
//   * more candidates than the output buffer holds  -> throws, with the count
//     the caller should have asked for (the DeviceAppendBuffer contract);
//   * a probe sequence exhausted, i.e. the table saturated -> throws.
//
// Neither can be "handled" by returning what was found: a truncated result from
// an exclusion-adjacent search is exactly the failure this project is built to
// avoid.
#pragma once

#include "core/collision.h"
#include "core/uint.h"

#include <cstdint>
#include <span>
#include <vector>

namespace pte::gpu {

struct ChampionParams {
  // Hash table size, as a power of two. The default is P4's measured shape:
  // 2^26 slots x 8 B = 512 MB, comfortably inside 8 GB alongside sq[] and the
  // staging buffers.
  int table_log = 26;

  // Occupancy the windows are sized to. Linear probing degrades sharply past
  // ~0.7; 0.5 is where P4 measured 1.0e9 inserts/s, and it sets the window cap
  // at slots * max_load = 3.35e7 pairs.
  double max_load = 0.5;

  // Capacity of the candidate output buffer, per window. Candidates are
  // collisions plus fingerprint false positives; at the default table size the
  // expected false-positive count is well under one per window (see the
  // arithmetic at kFingerprintBits in champion_impl.cu), so this is slack, not
  // a budget.
  std::size_t max_candidates = 1u << 16;
};

struct ChampionStats {
  std::uint64_t points = 0;
  std::uint64_t windows = 0;
  std::uint64_t pairs = 0;           // must come out at exactly C(D,2)
  std::uint64_t candidates = 0;      // what the device reported
  std::uint64_t false_positives = 0; // candidates the exact re-check rejected
  std::uint64_t hits = 0;

  // Where the wall clock went, split at the one seam that matters: the window
  // planner is an O(D) host sweep per window and the device sits idle through
  // it, so these two numbers say whether the search is device-bound or
  // planner-bound. They are wall time, not device time -- `device_seconds`
  // covers the H2D copies, the table clear, the launch and the synchronizing
  // read-back, which is what a window actually costs.
  double plan_seconds = 0.0;
  double device_seconds = 0.0;
};

// Every {a,b} | {c,d} with nu_a^2 + nu_b^2 == nu_c^2 + nu_d^2, as indices into
// `nu`, in canonical form and sorted -- the same answer cpu/champion's
// reference produces, by a different mechanism.
//
// `nu` must be the distinct values of a canonized point set (nu is injective
// there, which is what makes equal-key pairs automatically disjoint) and is
// limited to 2^20 points, because a pair of indices is packed into one 64-bit
// table slot alongside the fingerprint.
//
// Throws std::invalid_argument on a repeated nu, on more than 2^20 points, or
// on a malformed table size; std::runtime_error if a window overflows the
// candidate buffer or saturates the table; core::WindowIndivisible if one key
// value carries more pairs than a window can hold.
std::vector<core::PairCollision>
champion_meet(std::span<const core::U128> nu, const ChampionParams &params = {},
              ChampionStats *stats = nullptr);

} // namespace pte::gpu
