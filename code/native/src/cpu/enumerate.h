// Per-norm E_4 enumeration (search S-A of drafts/gpu_algorithms.md).
//
// The whole search reduces to one condition. For an admissible L1 -- a product
// of split primes with at least four representations as a sum of two squares --
// build the point set U(L1), take nu(z) = Re(z^2) for each canonized point, and
// look for four points splitting as {z1,z2} | {z3,z4} with
//
//     nu1^2 + nu2^2 == nu3^2 + nu4^2.
//
// That is the L_{2,b} form of the E_4 tree condition; see paper.md section 2
// and drafts/gpu_algorithms.md section 1.
#ifndef PTE_CPU_ENUMERATE_H_
#define PTE_CPU_ENUMERATE_H_

#include "gaussian.h"

#include <cstdint>
#include <functional>
#include <vector>

namespace pte::cpu {

// Above this bound on L1 the key nu1^2 + nu2^2 <= 2*L1^2 no longer fits in 64
// bits, and unsigned wraparound turns the meet from exact arithmetic into a
// filter (a sound one: equal keys always agree modulo 2^64, so there are no
// false negatives -- only false positives, which verify_exact rejects).
inline constexpr std::int64_t kExactKeyBound = 3037000499; // floor(2^63 ^ 1/2)

// One admissible norm: L1 together with the split primes multiplying to it,
// listed with multiplicity in ascending order.
struct Norm {
  std::int64_t l1 = 0;
  std::vector<std::uint64_t> primes;
};

// A candidate solution: four indices into a point set, as {a,b} | {c,d}.
struct Quadruple {
  int a = 0, b = 0, c = 0, d = 0;
};

// Split primes (p = 1 mod 4) up to and including `limit`, ascending.
std::vector<std::uint64_t> split_primes_up_to(std::uint64_t limit);

// Calls `cb` once per admissible L1 <= bound, in no particular order.
//
// Admissible means: a product of split primes only (an inert prime would divide
// every root, and the ramified 2 is removed by the odd-L1 normalization) whose
// point set has at least four elements, i.e. prod(b_i + 1) >= 8. Repeated
// primes are included -- roughly a fifth of admissible norms are non-squarefree
// and they can carry primitive solutions, so skipping them would silently break
// exhaustiveness.
void for_each_admissible_norm(std::int64_t bound,
                              const std::function<void(const Norm &)> &cb);

// Number of points a norm with these primes will have, without building them:
// ceil(prod(b_i + 1) / 2). Used to size buffers and to bucket work by D.
std::size_t point_count(const std::vector<std::uint64_t> &primes);

// The nu values of a norm's point set, in point-set order.
std::vector<std::int64_t> nu_values(const std::vector<Gauss> &points);

// Reference meet: every {a,b} | {c,d} with matching key, found by walking
// 4-subsets and testing the three ways to split each into two pairs. This is
// the CPU oracle the GPU kernel is checked against.
//
// No disjointness test is needed: nu is injective on a point set, so two
// distinct pairs with equal key cannot share a point.
std::vector<Quadruple> meet_reference(const std::vector<std::int64_t> &nu);

// Recomputes a candidate's key in 128-bit arithmetic. Above kExactKeyBound the
// device kernel works modulo 2^64, so every hit must pass through here.
bool verify_exact(const std::vector<std::int64_t> &nu, const Quadruple &q);

// The primitivity rule (drafts/gpu_algorithms.md A3.1). A quadruple whose four
// points are *all* divisible by some prime p is p^2 times a solution at the
// smaller norm L1/p^2 -- a duplicate that would otherwise be re-reported at
// every square multiple. A quadruple with only some points divisible by p is
// genuinely primitive and must be kept: that is exactly the case where one pair
// uses the rational factor p = q*conj(q) while another uses q^2.
//
// Only a repeated prime can divide a point at all: p | z needs both q and its
// conjugate to divide z, i.e. exponent b_p >= 2. Squarefree norms therefore
// always pass, and `primes` is the norm's factorization with multiplicity.
bool is_primitive(const std::vector<Gauss> &points,
                  const std::vector<std::uint64_t> &primes, const Quadruple &q);

} // namespace pte::cpu

#endif // PTE_CPU_ENUMERATE_H_
