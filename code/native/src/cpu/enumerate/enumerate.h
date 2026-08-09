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
#pragma once

#include "cpu/core/precision.h"
#include "cpu/core/solution.h"
#include "cpu/gaussian/gaussian.h"

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <span>
#include <vector>

namespace pte::cpu {

// An admissible norm together with its split-prime factorization -- not just a
// number, hence the name; `norm_t` is the bare integer width.
//
// The factorization is stored as *indices into the Enumerator's factor base*,
// ascending with multiplicity, rather than as prime values. The DFS knows the
// index at the moment it emits the norm, so storing values would mean throwing
// that away and searching the base again to recover it. Enumerator::prime()
// and Enumerator::primes_of() resolve indices when a caller needs values.
struct DecomposedNorm {
  norm_t l1 = 0;
  std::vector<factor_index_t> factors;
};

// Split primes (p = 1 mod 4) up to and including `limit`, ascending.
std::vector<factor_t> split_primes_up_to(factor_t limit);

// Number of points a norm with this factorization will have, without building
// them: ceil(prod(b_i + 1) / 2). Used to size buffers and to bucket work by D.
//
// Templated on the range because it only needs to spot runs of equal entries,
// which works identically on prime values and on factor-base indices.
template <typename Range> std::size_t point_count(const Range &factors) {
  std::size_t divisors = 1;
  std::size_t run = 0;
  const auto n = static_cast<std::size_t>(std::size(factors));
  for (std::size_t i = 0; i < n; ++i) {
    ++run;
    if (i + 1 == n || std::data(factors)[i + 1] != std::data(factors)[i]) {
      divisors *= run + 1;
      run = 0;
    }
  }
  return (divisors + 1) / 2;
}

// Drives the enumeration for one bound, owning everything that is shared across
// norms.
//
// This exists because the free-function form recomputed too much. Every norm
// re-derived the Cornacchia decomposition of each of its prime factors, which
// measured at **60% of point-set construction time** and a 41x-62x redundancy
// factor (36,675 calls against 891 distinct primes at B = 10^6; 436,126 against
// 7,072 at B = 10^7). The factor base and its decompositions are properties of
// the bound, not of a norm, so they belong to an object that outlives one.
//
// **This is the reference enumerator, not the production one.** It rebuilds
// every point set from {1,0} -- 137 ns/norm against the Spine's 54 -- and it
// walks the tree smallest-prime-first, which puts 89% of the work in one
// segment and cannot be parallelized. cpu/generate/ is what a campaign runs.
//
// It stays because it is obviously correct and independently written: the Spine
// is checked against it norm for norm, which is worth more than the speed it
// gives up. Do not refactor the two to share code.
class Enumerator {
public:
  // Builds the factor base and decomposes every prime in it exactly once.
  // Throws std::invalid_argument on a negative bound, and NotASplitPrime if the
  // sieve ever yields a non-split prime, which would mean the sieve is wrong.
  //
  // A bound below the smallest admissible norm (5*13*17 = 1105) is *not* an
  // error: "no admissible norms exist below 1000" is a true and useful answer,
  // so the enumeration is simply empty.
  explicit Enumerator(norm_t bound);

  norm_t bound() const { return bound_; }

  // The split primes that can divide an admissible norm <= bound(), ascending.
  const std::vector<factor_t> &factor_base() const { return primes_; }

  // Resolve a factorization index to its prime, and its decomposition.
  factor_t prime(factor_index_t i) const { return primes_[i]; }
  const Gauss &prime_point(factor_index_t i) const { return prime_points_[i]; }

  // The prime values of a norm's factorization, with multiplicity. Only needed
  // by callers that must do arithmetic with the primes themselves, such as
  // is_primitive; the enumeration itself stays on indices.
  std::vector<factor_t> primes_of(const DecomposedNorm &norm) const;

  // Calls `cb(const DecomposedNorm&)` once per admissible L1 <= bound(), in no
  // particular order.
  //
  // Admissible means: a product of split primes only (an inert prime would
  // divide every root, and the ramified 2 is removed by the odd-L1
  // normalization) whose point set has at least four elements, i.e.
  // prod(b_i + 1) >= 8. Repeated primes are included -- roughly a fifth of
  // admissible norms are non-squarefree and they can carry primitive solutions,
  // so skipping them would silently break exhaustiveness.
  //
  // The callback is a template parameter rather than a std::function so the DFS
  // body can inline it; this loop runs ~2e10 times over a full campaign.
  //
  // The DecomposedNorm handed to `cb` aliases the DFS's working state and is
  // invalidated as soon as `cb` returns. Copy it if you intend to keep it.
  template <typename F> void for_each_norm(F &&cb) const;

  // The point set of an admissible norm, using the memoized decompositions.
  // A direct indexed lookup -- no search, which is the point of storing
  // indices in DecomposedNorm rather than values.
  std::vector<Gauss> point_set(const DecomposedNorm &norm) const;

private:
  norm_t bound_ = 0;
  std::vector<factor_t> primes_;
  std::vector<Gauss> prime_points_; // primes_[i] = N(prime_points_[i])
};

// Map a point set through `f`, preserving order. `nu_values` is the one use
// that matters today; the generic form exists because the meet will want norms
// and 2-adic valuations through the same loop.
template <typename F>
auto transform_points(const std::vector<Gauss> &points, F &&f) {
  std::vector<decltype(f(points.front()))> out;
  out.reserve(points.size());
  for (const Gauss &z : points)
    out.push_back(f(z));
  return out;
}

// The nu values of a norm's point set, in point-set order.
inline std::vector<norm_t> nu_values(const std::vector<Gauss> &points) {
  return transform_points(
      points, [](const Gauss &z) { return static_cast<norm_t>(z.nu()); });
}

// The meet *oracle*: every {a,b} | {c,d} with matching key, by walking every
// 4-subset and testing all three ways to split it into two pairs.
//
// "Oracle", not "the CPU implementation", and the distinction is the whole
// reason this function exists:
//
//   - It is exact throughout (exact_key_t), never a filter. The production
//     meet compares modulo 2^64 and needs confirm_candidate behind it; this
//     does not, so it can adjudicate what the fast paths claim.
//   - It is deliberately the naive O(D^4) form, with no bucketing, no
//     precomputed squares and no early exit -- so it shares no structure with
//     the optimized paths and cannot share a bug with them.
//   - It is written to be read, not to be fast. At 2.0 ns per partition it is
//     usable to B ~ 10^7 and no further.
//
// The production CPU meet, when it arrives, will be the D-bucketed
// filter-plus-confirm path and will be checked *against* this. Keeping the two
// separate is what makes the cross-check evidence rather than a tautology.
//
// No disjointness test is needed: nu is injective on a canonized point set, so
// two distinct pairs with equal key cannot share a point.
std::vector<E4Sol> meet_oracle(const std::vector<norm_t> &nu);

// Confirms a candidate by recomputing its key in exact arithmetic.
//
// **This is a production step, not a test helper.** Above kExactKeyBound the
// device kernel compares keys modulo 2^64 and therefore emits candidates; this
// is the layer-2 arbiter from precision.h that turns a candidate into a result.
// Skipping it does not make the search faster in any way that matters -- there
// are a handful of candidates -- it just makes the output wrong, silently, and
// only above L1 = 3.037e9 where nothing else would notice.
//
// Also rejects malformed quadruples (out-of-range or repeated indices), so a
// garbled hit from a buggy kernel cannot be mistaken for a solution.
bool confirm_candidate(std::span<const norm_t> nu, const E4Sol &q);

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
bool is_primitive(std::span<const Gauss> points, std::span<const factor_t> primes,
                  const E4Sol &q);

// ---------------------------------------------------------------------------
// Implementation of the templated member
// ---------------------------------------------------------------------------

template <typename F> void Enumerator::for_each_norm(F &&cb) const {
  if (primes_.empty())
    return;

  DecomposedNorm norm;
  norm.factors.reserve(16);

  // Walks primes in ascending index order, choosing an exponent for each, so
  // every multiset of split primes with product <= bound is visited once.
  const auto dfs = [&](auto &&self, std::size_t start, norm_t value,
                       std::size_t divisors) -> void {
    for (std::size_t i = start; i < primes_.size(); ++i) {
      const auto p = static_cast<norm_t>(primes_[i]);
      if (value > bound_ / p)
        break; // primes ascend, so nothing further fits either

      norm_t v = value * p;
      std::size_t exponent = 1;
      while (true) {
        norm.factors.push_back(static_cast<factor_index_t>(i));
        const std::size_t d = divisors * (exponent + 1);
        if (d >= 8) { // >= 4 points, the minimum for an E_4
          norm.l1 = v;
          cb(static_cast<const DecomposedNorm &>(norm));
        }
        self(self, i + 1, v, d);

        if (v > bound_ / p)
          break;
        v *= p;
        ++exponent;
      }
      norm.factors.resize(norm.factors.size() - exponent);
    }
  };
  dfs(dfs, 0, 1, 1);
}

} // namespace pte::cpu
