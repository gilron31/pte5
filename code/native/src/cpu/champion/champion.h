// S-E, the champion meet: collisions inside one enormous point set.
//
// ---------------------------------------------------------------------------
// What search this is, and what it is for
// ---------------------------------------------------------------------------
//
// The exhaustive campaign (S-A) sweeps every admissible L1 below a bound and
// answers an *exclusion* question. S-E does the opposite: it takes one
// hand-chosen norm with an enormous point set -- 20 distinct split primes give
// D = 2^19 at L1 ~ 2.5e36, far outside any exhaustive bound -- and asks whether
// that one norm contains a collision. It supports **discovery, not exclusion**,
// and it is the search drafts/gpu_implementation_plan.md ranks first for the
// device.
//
// The question is the same one as ever (gpu_algorithms.md 1):
//
//     nu_a^2 + nu_b^2 == nu_c^2 + nu_d^2   over disjoint pairs {a,b}, {c,d}
//
// Only the scale changes, and the scale changes everything about the
// implementation:
//
//   D          2^19 rather than <= 192, so C(D,2) = 1.374e11 keys
//   nu         122 bits rather than 40, so the key is ~244 and needs core::UInt
//   the meet   neither meet_direct (D^4) nor meet_hashed (all keys resident)
//              can run: the key array alone would be 2.2 TB
//
// So the algorithm is: partition the *key range* into windows small enough to
// hold (core/window.h), and inside each window find equal keys. That is the
// 1000x correction of gpu_implementation_plan.md 3.1 -- windowing by key value
// rather than by hash means each pair is enumerated once instead of R times.
//
// ---------------------------------------------------------------------------
// This file is the reference; the device kernel is the fast one
// ---------------------------------------------------------------------------
//
// gpu/champion/champion.h answers the same question with an atomicCAS hash
// table, at ~1e9 inserts/s. This one sorts each window's keys and reads off
// equal-key runs -- slower by a wide margin, and **exact by construction**: the
// keys it compares are full 256-bit values, so it has no false positives to
// confirm and no fingerprint to collide. It is the arbiter of the two, in the
// same relation as meet_oracle to the production meet: written to be obviously
// right rather than fast, deliberately sharing no collision-detection code with
// the implementation it checks.
//
// Both use the same core/window.h partition. That is on purpose and is the one
// piece they share: value-windowing is not a heuristic, it is the definition of
// which pairs exist, and having two versions of it would mean the cross-check
// could not tell "we disagree" from "we enumerated different pairs". It is
// tested on its own terms in core/window_test.cc against brute force.
//
// ---------------------------------------------------------------------------
// Precondition, and why it is checked
// ---------------------------------------------------------------------------
//
// nu must be **distinct** -- the values of a canonized point set, where nu is
// injective (gaussian.h). Everything downstream rests on it: injectivity is
// what makes two pairs with equal key automatically disjoint, which is what
// licenses the absence of a disjointness test here and in every kernel. A
// caller who passes raw products rather than canonized points violates it, so
// champion_meet checks and throws rather than returning quadruples that share a
// point.
#pragma once

#include "core/collision.h"
#include "core/uint.h"

#include <cstdint>
#include <span>
#include <vector>

namespace pte::cpu {

struct ChampionParams {
  // Cap on pairs per value window. This is a *memory* knob: the reference holds
  // one entry (key + two indices, 40 B) per pair in flight, so the default is
  // 40 MB. The device path sizes it from its hash table instead.
  std::uint64_t max_pairs_per_window = 1u << 20;
};

// Counters worth having on a run that takes minutes, and worth asserting on in
// a test: `pairs` must come out at exactly C(D,2), which is the statement that
// no pair was skipped.
struct ChampionStats {
  std::uint64_t points = 0;
  std::uint64_t windows = 0;
  std::uint64_t pairs = 0;
  std::uint64_t hits = 0;
};

// Every {a,b} | {c,d} with nu_a^2 + nu_b^2 == nu_c^2 + nu_d^2, as indices into
// `nu`, in canonical form and sorted.
//
// The result is exact: no confirmation step is needed or offered, because
// nothing here compares a truncated key.
//
// Throws std::invalid_argument if `nu` contains a repeated value (see the
// precondition above), and core::WindowIndivisible if a single key value
// carries more pairs than max_pairs_per_window allows.
std::vector<core::PairCollision>
champion_meet(std::span<const core::U128> nu, const ChampionParams &params = {},
              ChampionStats *stats = nullptr);

} // namespace pte::cpu
