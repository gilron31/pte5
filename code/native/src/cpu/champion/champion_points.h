// Point sets for the champion regime, where the coordinates no longer fit the
// enumeration's types.
//
// ---------------------------------------------------------------------------
// Why this is not gaussian.h's point_set
// ---------------------------------------------------------------------------
//
// cpu/gaussian's point_set is built on GaussT<coord_t> with coord_t = int32,
// which is exactly right at the design bound: L1 <= 10^12 gives |m|,|n| <= 10^6
// and 8 bytes per point (cpu/core/precision.h). A champion norm is
// L1 ~ 2.5e36, so the coordinates are 61 bits and every intermediate product is
// 122 -- the int32 path does not merely overflow, it overflows silently, which
// is the failure precision.h exists to prevent.
//
// So the multiplicative construction is repeated here over mpz_class. That is
// deliberately a *different* implementation rather than a template
// instantiation of the same one: GaussT's `wide` result type has no meaning at
// arbitrary precision, and there is nothing to be gained by making one
// construction serve widths a hundred thousand times apart. The check that they
// agree is a test, not a shared header -- champion_points_test.cc runs both on
// small norms and requires identical nu values.
//
// The cost is one-off and small against the meet it feeds: a 20-prime champion
// is ~2^20 GMP multiplications of 122-bit values to produce D = 2^19 points,
// which is a fraction of a second against S-E's minutes.
//
// ---------------------------------------------------------------------------
// What "champion" means for the prime set
// ---------------------------------------------------------------------------
//
// D = ceil(prod(b_i + 1) / 2), so k *distinct* split primes give D = 2^(k-1),
// and D is what the search is quadratic in. Twenty distinct split primes is
// D = 2^19 at L1 = 2.5e36 -- open question 5 of the plan is which twenty, and
// `smallest_split_primes` answers only the obvious version of it (take the
// smallest, which maximises D per unit of L1).
#pragma once

#include "core/uint.h"
#include "cpu/core/precision.h"

#include <cstdint>
#include <gmpxx.h>
#include <vector>

namespace pte::cpu {

// A champion norm: its prime factorization, the norm itself, and the point set
// as the nu values the meet consumes.
struct ChampionNorm {
  std::vector<factor_t> primes; // distinct split primes, ascending
  mpz_class l1;                 // their product
  std::vector<core::U128> nu;   // nu(z) = m^2 - n^2 over canonized points
};

// The `count` smallest primes congruent to 1 mod 4, which are the split primes
// (gaussian.h). Ascending.
//
// Smallest-first maximises D for a given L1, which is the only sizing rule the
// plan commits to; whether a champion should instead be tuned to the section-4
// ladder congruences is open question 5 and is not decided here.
std::vector<factor_t> smallest_split_primes(int count);

// Builds U(L1) for the product of `primes`, canonized to the first octant and
// deduplicated, and returns it as the champion norm above.
//
// `primes` must be distinct split primes (p = 1 mod 4). Repeated primes are
// rejected rather than handled: the exponent-sweeping construction they need is
// the campaign's business (cpu/gaussian), the champion regime has no reason to
// want them, and accepting them silently would make D disagree with 2^(k-1).
//
// Throws std::invalid_argument on a repeat or on more than 21 primes (D = 2^20
// is the ceiling gpu/champion can index), and NotASplitPrime on a prime that is
// not 1 mod 4.
//
// The nu values are distinct by construction -- nu is injective on canonized
// points -- so the result is directly the precondition champion_meet requires.
ChampionNorm champion_norm(const std::vector<factor_t> &primes);

} // namespace pte::cpu
