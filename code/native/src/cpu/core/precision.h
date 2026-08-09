// The precision policy for the enumeration.
//
// This project's output is an *exclusion* claim -- "no E_5 tower exists below
// this bound" -- and the characteristic way such a claim dies is silent
// overflow: a search that quietly stops finding things and reports nothing,
// which is indistinguishable from a search that found nothing. So arithmetic
// width is a first-class design decision here rather than an implementation
// detail, and it is stated once, in this file, instead of being rediscovered as
// a `using u64` at the top of each .cc.
//
// ---------------------------------------------------------------------------
// Three layers, three widths
// ---------------------------------------------------------------------------
//
//   1. ENUMERATION -- fixed width, and allowed to be a *filter*.
//      Point sets, nu values, and the GPU meet. May wrap, on one condition:
//      wrapping must produce false positives only, never false negatives. The
//      meet satisfies this because equal keys are always congruent mod 2^64.
//
//   2. CONFIRMATION -- exact, and the arbiter.
//      Every candidate the filter emits is re-checked here in `exact_key_t`
//      before it is called a solution. This is the layer that makes layer 1's
//      sloppiness sound, so it must never be skipped or compiled out.
//
//   3. CERTIFICATE -- arbitrary precision, GMP.
//      Checking a recovered PTE identity or a full E_5 tower. Not yet written;
//      see the bit count below for why it cannot reuse layer 2.
//
// ---------------------------------------------------------------------------
// The numbers, at the design bound B = kMaxNorm = 10^12
// ---------------------------------------------------------------------------
//
//   quantity                       bound     bits   type
//   ---------------------------------------------------------------
//   Gauss coordinates m, n         sqrt(B)     20   coord_t
//   L1, nu = m^2 - n^2             B           40   norm_t
//   pair key nu_a^2 + nu_b^2       2B^2        81   exact_key_t
//   PTE power sums, k <= 15        16*B^7.5   311   GMP  <-- layer 3
//
// The 81-bit pair key is the load-bearing row: it does not fit in 64 bits, but
// it fits in __int128 with 46 bits to spare. Concretely, `exact_key_t` stays
// exact for L1 < 9.2e18, which is nine million times past B. **GMP is not
// needed anywhere in layers 1 or 2**, and putting it there would cost a heap
// allocation in the one loop that must not have one.
//
// GMP becomes unavoidable only at layer 3: verifying an ideal symmetric PTE of
// size 16 means power sums to k = 15 of sixteen integers of size ~sqrt(2*L1),
// which is ~311 bits. A size-32 identity would need ~638.
#pragma once

#include <cstdint>
#include <limits>

namespace pte::cpu {

// The design bound on L1. Everything below is sized against it; raising it is
// a deliberate act that the static_asserts at the bottom will police.
inline constexpr std::int64_t kMaxNorm = 1'000'000'000'000; // 10^12

// ---------------------------------------------------------------------------
// Layer 1: enumeration
// ---------------------------------------------------------------------------

// Gauss coordinates. int64 rather than the int32 that would suffice (20 bits at
// B = 10^12), because nothing currently transfers raw coordinates -- the meet
// moves nu values, which need 40 bits. Narrowing GaussT<std::int32_t> halves
// the bytes per point and becomes worth doing when device-resident generation
// puts coordinates on the bus; the template is already shaped for it.
using coord_t = std::int64_t;

// L1 and nu(z). Both bounded by kMaxNorm.
using norm_t = std::int64_t;

// ---------------------------------------------------------------------------
// Layer 2: confirmation
// ---------------------------------------------------------------------------

// The exact pair key. __int128 is a GNU extension; the opt-in keeps -Wpedantic
// quiet, and GCC/Clang/nvcc is the entire universe this project targets.
__extension__ using exact_key_t = __int128;

// Above this bound on L1, the key nu_a^2 + nu_b^2 <= 2*L1^2 no longer fits in
// an unsigned 64-bit word, and the GPU meet's wraparound turns from exact
// arithmetic into a filter -- a sound one (equal keys stay congruent, so there
// are no false negatives), but one whose false positives layer 2 must reject.
//
// = floor(sqrt(2^63)), i.e. the largest L1 with 2*L1^2 < 2^64.
inline constexpr norm_t kExactKeyBound = 3'037'000'499;

// True if a batch of norms up to `max_l1` keeps the meet in exact arithmetic,
// so a caller can tell the two regimes apart rather than crossing the boundary
// unawares. Below it, a candidate is already a solution; above it, confirmation
// is doing real work.
constexpr bool meet_is_exact(norm_t max_l1) { return max_l1 <= kExactKeyBound; }

// ---------------------------------------------------------------------------
// The asserts that make the table above a contract rather than a comment
// ---------------------------------------------------------------------------

static_assert(sizeof(exact_key_t) == 16,
              "the pair key needs 81 bits at B = 10^12; 64 is not enough");

// A coordinate must hold sqrt(kMaxNorm)...
static_assert(static_cast<double>(std::numeric_limits<coord_t>::max()) >
                  1.1e6, // sqrt(10^12), with headroom
              "coord_t cannot hold sqrt(kMaxNorm)");

// ...a norm must hold kMaxNorm itself...
static_assert(std::numeric_limits<norm_t>::max() >= kMaxNorm,
              "norm_t cannot hold kMaxNorm");

// ...and the exact key must hold 2*kMaxNorm^2 with room to spare.
static_assert(static_cast<exact_key_t>(kMaxNorm) * kMaxNorm * 2 > 0,
              "exact_key_t overflows on 2*kMaxNorm^2");

// The filter boundary is genuinely inside the exact type's range: if this ever
// failed, the "confirm in wider arithmetic" story would be circular.
static_assert(static_cast<exact_key_t>(kExactKeyBound) * kExactKeyBound * 2 > 0,
              "exact_key_t cannot confirm candidates at kExactKeyBound");

} // namespace pte::cpu
