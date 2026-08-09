// Gaussian-integer arithmetic for the E_4 enumeration.
//
// The enumeration works at the "norm-L1 level": a root pair of an E_4 solution
// is recovered from a point z = m + ni with N(z) = L1 by the twist
// z -> (1+i)z, which sends norm L1 to norm 2L1 and gives the root pair
// (m-n, m+n). Staying at norm L1 keeps every operand one bit smaller, so the
// twist is applied once, at output time.
//
// See drafts/gpu_algorithms.md section 1 for the derivation, and
// cpu/core/precision.h for why the widths here are what they are.
#pragma once

#include "cpu/core/precision.h"

#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace pte::cpu {

// A Gaussian integer, templated on the coordinate type.
//
// The template parameter is not decoration: at B = 10^12 a coordinate needs 20
// bits, so GaussT<std::int32_t> is valid and halves the bytes per point. That
// matters once point coordinates are what crosses the PCIe bus (device-resident
// generation); until then the concrete alias below stays 64-bit, which cannot
// overflow for any bound we can enumerate.
//
// The subtlety worth stating, because getting it wrong is exactly the silent
// overflow precision.h exists to prevent: **norm() and nu() do not fit in T.**
// They reach L1 = 10^12 while a coordinate reaches 10^6. Hence `wide`, and
// hence the products below being formed in `wide` before any narrowing.
template <typename T> struct GaussT {
  static_assert(std::is_signed_v<T> && std::is_integral_v<T>,
                "GaussT needs a signed integral coordinate type");

  // Wide enough for a product of two coordinates: 2*bits(T), so the norm of a
  // GaussT<int32_t> lands in int64 and the norm of a GaussT<int64_t> in
  // __int128.
  using wide = std::conditional_t<sizeof(T) <= 4, std::int64_t, exact_key_t>;

  T re = 0;
  T im = 0;

  friend bool operator==(const GaussT &a, const GaussT &b) {
    return a.re == b.re && a.im == b.im;
  }
  friend bool operator!=(const GaussT &a, const GaussT &b) { return !(a == b); }
  // Lexicographic; only needed to give point sets a deterministic order.
  friend bool operator<(const GaussT &a, const GaussT &b) {
    return a.re != b.re ? a.re < b.re : a.im < b.im;
  }

  GaussT conj() const { return {re, static_cast<T>(-im)}; }

  // N(z) = m^2 + n^2, and nu(z) = Re(z^2) = m^2 - n^2.
  //
  // nu is *injective* on a point set: m^2 + n^2 = L1 together with m^2 - n^2 =
  // nu determines (m,n). That is why a key collision between two distinct pairs
  // is automatically disjoint, which removes the disjointness test from the
  // meet entirely -- on both the CPU and the GPU.
  wide norm() const { return wide{re} * re + wide{im} * im; }
  wide nu() const { return wide{re} * re - wide{im} * im; }

  // Representative of z under the action of the units and conjugation: the
  // point of its orbit lying in the first octant, m >= n >= 0. The E_4
  // invariants are all constant on these orbits, so canonizing is what makes a
  // "point set" a set.
  GaussT canonize() const {
    const T a = re < 0 ? static_cast<T>(-re) : re;
    const T b = im < 0 ? static_cast<T>(-im) : im;
    return a >= b ? GaussT{a, b} : GaussT{b, a};
  }

  // True if the rational prime p divides z, i.e. both coordinates. Only a
  // repeated prime can divide a point at all; see is_primitive.
  bool divisible_by(T p) const { return re % p == 0 && im % p == 0; }

  // Products are formed in `wide` and narrowed back. The narrowing is safe for
  // any product of two points whose norms multiply to at most kMaxNorm, which
  // is every multiplication the enumeration performs: |a.re * b.re| <=
  // sqrt(N(a) * N(b)) = sqrt(N(ab)).
  friend GaussT operator*(const GaussT &a, const GaussT &b) {
    return {static_cast<T>(wide{a.re} * b.re - wide{a.im} * b.im),
            static_cast<T>(wide{a.re} * b.im + wide{a.im} * b.re)};
  }

  // Both z*q and z*conj(q) from four integer multiplies instead of eight.
  //
  // Deliberately not an operator and deliberately not two calls to operator*:
  // it returns two results from shared subexpressions, which is the whole point
  // of it. This is the inner operation of point-set construction.
  // GILRO: use lvalue ref not ptrs.
  void mul_and_conj_mul(const GaussT &q, GaussT *out, GaussT *out_conj) const {
    const wide x = wide{re} * q.re;
    const wide y = wide{re} * q.im;
    const wide w = wide{im} * q.re;
    const wide v = wide{im} * q.im;
    *out = {static_cast<T>(x - v), static_cast<T>(y + w)};      // z * q
    *out_conj = {static_cast<T>(x + v), static_cast<T>(w - y)}; // z * conj(q)
  }
};

// The concrete point type. See coord_t in precision.h for the width choice.
using Gauss = GaussT<coord_t>;

// Floor of the square root of n, exact for every 64-bit input.
//
// Exposed rather than kept private to the .cc because it is the one piece of
// hand-rolled numerics in this layer, it is used inside Cornacchia's descent
// where a wrong answer would silently drop a prime, and an earlier version of
// it in this project was in fact wrong. gaussian_test.cc checks it against
// GMP's mpz_sqrt over random and boundary inputs.
//
// Kept rather than simply calling mpz_sqrt: this runs once per prime, and an
// mpz init plus library call costs several times more than the one or two
// Newton corrections here.
std::uint64_t isqrt(std::uint64_t n);

// True if n is prime, by deterministic Miller-Rabin. Exact over the whole
// 64-bit range, not probabilistic. Exposed for the same reason as isqrt.
bool is_prime(std::uint64_t n);

// Thrown when a caller hands decompose_prime something that is not a prime
// congruent to 1 mod 4. A distinct type because the enumeration's factor base
// contains only split primes by construction, so this is a bug at the call
// site, not an expected outcome -- and swallowing it would silently delete
// every solution at every norm divisible by that prime.
class NotASplitPrime : public std::invalid_argument {
public:
  using std::invalid_argument::invalid_argument;
};

// Writes p = a^2 + b^2 with a >= b > 0 for a prime p = 1 (mod 4), by
// Tonelli-Shanks for sqrt(-1) followed by Cornacchia's descent.
//
// Throws NotASplitPrime for p = 2 (ramified), for p = 3 (mod 4) (inert), and
// for composites. It does not return a status: an exclusion search must never
// be able to treat "I could not factor this" as "there is nothing here".
//
// Primality is checked by deterministic Miller-Rabin, not inferred from the
// descent succeeding. A composite that is itself a product of split primes
// (85 = 5*17, 221 = 13*17) *is* a sum of two squares, so the descent closes and
// hands back a valid decomposition of a non-prime -- and those are precisely
// the composites a caller is likely to pass by accident.
Gauss decompose_prime(std::uint64_t p);

// All Gaussian integers of norm `prod(primes)`, canonized and deduplicated,
// where `primes` lists split primes with multiplicity. The result has
// ceil(prod(b_i + 1) / 2) elements for L1 = prod q_i^{b_i}.
//
// Built multiplicatively: branching over *each copy* of a repeated prime is
// what realizes every exponent 0 <= c <= b_i, which is the non-squarefree
// coverage that any exhaustiveness claim rests on (paper.md Prop. 5.1).
//
// Recomputes every prime's decomposition on every call, which measured at 60%
// of point-set construction. Enumerator::point_set memoizes and is what the
// enumeration should use; this overload remains for tests and one-off queries.
std::vector<Gauss> point_set(const std::vector<std::uint64_t> &primes);

// As above, but taking decompositions the caller already has, one per entry of
// `primes`. This is the form the enumeration uses.
std::vector<Gauss> point_set_from(const std::vector<Gauss> &prime_points);

} // namespace pte::cpu
