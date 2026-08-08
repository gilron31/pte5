// Gaussian-integer arithmetic for the E_4 enumeration.
//
// The enumeration works at the "norm-L1 level": a root pair of an E_4 solution
// is recovered from a point z = m + ni with N(z) = L1 by the twist
// z -> (1+i)z, which sends norm L1 to norm 2L1 and gives the root pair
// (m-n, m+n). Staying at norm L1 keeps every operand one bit smaller, so the
// twist is applied once, at output time.
//
// See drafts/gpu_algorithms.md section 1 for the derivation.
#ifndef PTE_CPU_GAUSSIAN_H_
#define PTE_CPU_GAUSSIAN_H_

#include <cstdint>
#include <vector>

namespace pte::cpu {

// A Gaussian integer. 64-bit coordinates are ample: an admissible L1 <= 1e12
// has |m|, |n| < 1e6.
struct Gauss {
  std::int64_t re = 0;
  std::int64_t im = 0;

  friend bool operator==(const Gauss &a, const Gauss &b) {
    return a.re == b.re && a.im == b.im;
  }
  friend bool operator!=(const Gauss &a, const Gauss &b) { return !(a == b); }
  // Lexicographic; only needed to give point sets a deterministic order.
  friend bool operator<(const Gauss &a, const Gauss &b) {
    return a.re != b.re ? a.re < b.re : a.im < b.im;
  }
};

inline Gauss conj(const Gauss &z) { return {z.re, -z.im}; }

inline Gauss mul(const Gauss &a, const Gauss &b) {
  return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}

inline std::int64_t norm(const Gauss &z) { return z.re * z.re + z.im * z.im; }

// Both z*q and z*conj(q) from four integer multiplies instead of eight. This is
// the inner operation of point-set construction, so the saving is worth the
// slightly awkward signature.
inline void mul_and_conj_mul(const Gauss &z, const Gauss &q, Gauss *out,
                             Gauss *out_conj) {
  const std::int64_t x = z.re * q.re;
  const std::int64_t y = z.re * q.im;
  const std::int64_t w = z.im * q.re;
  const std::int64_t v = z.im * q.im;
  *out = {x - v, y + w};      // z * q
  *out_conj = {x + v, w - y}; // z * conj(q)
}

// Representative of z under the action of the units and conjugation: the point
// of its orbit lying in the first octant, m >= n >= 0. The E_4 invariants are
// all constant on these orbits, so canonizing is what makes a "point set" a
// set.
inline Gauss canonize(const Gauss &z) {
  const std::int64_t a = z.re < 0 ? -z.re : z.re;
  const std::int64_t b = z.im < 0 ? -z.im : z.im;
  return a >= b ? Gauss{a, b} : Gauss{b, a};
}

// nu(z) = Re(z^2) = m^2 - n^2. On canonized points this is non-negative, and it
// is *injective*: m^2 + n^2 = L1 together with m^2 - n^2 = nu determines (m,n).
// Injectivity is why a key collision between two distinct pairs is
// automatically disjoint, which removes the disjointness test from the meet.
inline std::int64_t nu(const Gauss &z) { return z.re * z.re - z.im * z.im; }

// Writes p = a^2 + b^2 with a >= b > 0 for a prime p = 1 (mod 4), by
// Tonelli-Shanks for sqrt(-1) followed by Cornacchia's descent. Returns false
// if p is not such a prime (including p = 2 and the inert p = 3 (mod 4)).
bool decompose_prime(std::uint64_t p, Gauss *out);

// All Gaussian integers of norm `prod(primes)`, canonized and deduplicated,
// where `primes` lists split primes with multiplicity. The result has
// ceil(prod(b_i + 1) / 2) elements for L1 = prod q_i^{b_i}.
//
// Built multiplicatively: branching over *each copy* of a repeated prime is
// what realizes every exponent 0 <= c <= b_i, which is the non-squarefree
// coverage that any exhaustiveness claim rests on (paper.md Prop. 5.1).
std::vector<Gauss> point_set(const std::vector<std::uint64_t> &primes);

} // namespace pte::cpu

#endif // PTE_CPU_GAUSSIAN_H_
