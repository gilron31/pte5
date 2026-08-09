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
#include "cpu/gaussian/gaussian.h"

#include <cstdint>
#include <vector>

namespace pte::cpu {

// One admissible norm: L1 together with the split primes multiplying to it,
// listed with multiplicity in ascending order.
struct Norm {
  norm_t l1 = 0;
  std::vector<std::uint64_t> primes;
};

// A candidate solution: four indices into a point set, as {a,b} | {c,d}.
struct Quadruple {
  int a = 0, b = 0, c = 0, d = 0;
};

// Split primes (p = 1 mod 4) up to and including `limit`, ascending.
std::vector<std::uint64_t> split_primes_up_to(std::uint64_t limit);

// Number of points a norm with these primes will have, without building them:
// ceil(prod(b_i + 1) / 2). Used to size buffers and to bucket work by D.
std::size_t point_count(const std::vector<std::uint64_t> &primes);

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
// Not yet done, and the reason this is a class rather than a memo table: a DFS
// node's point set is its parent's doubled by one prime, so the point sets
// themselves can be carried down the spine instead of rebuilt per norm. That is
// A2 (largest-prime-outer) in drafts/gpu_algorithms.md, and this is where it
// will live.
class Enumerator {
public:
  // Builds the factor base and decomposes every prime in it exactly once.
  // Throws NotASplitPrime if the sieve ever yields a non-split prime, which
  // would mean the sieve itself is wrong.
  explicit Enumerator(norm_t bound);

  norm_t bound() const { return bound_; }

  // The split primes that can divide an admissible norm <= bound(), ascending.
  const std::vector<std::uint64_t> &factor_base() const { return primes_; }

  // Calls `cb(const Norm&)` once per admissible L1 <= bound(), in no particular
  // order.
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
  // The Norm handed to `cb` aliases the DFS's working state and is invalidated
  // as soon as `cb` returns. Copy it if you intend to keep it.
  template <typename F> void for_each_norm(F &&cb) const;

  // The point set of an admissible norm, using the memoized decompositions.
  std::vector<Gauss> point_set(const Norm &norm) const;

private:
  // Index of `p` in primes_, for the memo lookup. The DFS knows the index
  // already; this is for callers who arrive with a factorization.
  std::size_t index_of(std::uint64_t p) const;

  norm_t bound_ = 0;
  std::vector<std::uint64_t> primes_;
  std::vector<Gauss> prime_points_; // primes_[i] = N(prime_points_[i])
};

// The nu values of a norm's point set, in point-set order.
std::vector<norm_t> nu_values(const std::vector<Gauss> &points);

// Reference meet: every {a,b} | {c,d} with matching key, found by walking
// 4-subsets and testing the three ways to split each into two pairs. Exact
// throughout (exact_key_t), so this is the oracle the GPU kernel is checked
// against rather than a filter.
//
// No disjointness test is needed: nu is injective on a point set, so two
// distinct pairs with equal key cannot share a point.
std::vector<Quadruple> meet_reference(const std::vector<norm_t> &nu);

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
bool confirm_candidate(const std::vector<norm_t> &nu, const Quadruple &q);

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

// ---------------------------------------------------------------------------
// Implementation of the templated member
// ---------------------------------------------------------------------------

template <typename F> void Enumerator::for_each_norm(F &&cb) const {
  if (primes_.empty())
    return;

  Norm norm;
  norm.primes.reserve(16);

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
        norm.primes.push_back(primes_[i]);
        const std::size_t d = divisors * (exponent + 1);
        if (d >= 8) { // >= 4 points, the minimum for an E_4
          norm.l1 = v;
          cb(static_cast<const Norm &>(norm));
        }
        self(self, i + 1, v, d);

        if (v > bound_ / p)
          break;
        v *= p;
        ++exponent;
      }
      norm.primes.resize(norm.primes.size() - exponent);
    }
  };
  dfs(dfs, 0, 1, 1);
}

} // namespace pte::cpu
