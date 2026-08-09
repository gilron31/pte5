// The split-prime factor base: every prime p = 1 (mod 4) up to P_max, together
// with its Gaussian decomposition p = a^2 + b^2.
//
// This is one-time work (drafts/implementation_design.md 3.1). Measured
// single-threaded, segmented:
//
//   P_max      split primes    sieve    Cornacchia    total    resident
//   1.54e8       4,327,911     0.73 s      3.86 s     4.6 s     0.1 MB
//   1.54e9      38,308,127     7.33 s     38.74 s    46.1 s     0.2 MB
//   1.54e10    343,618,258    74.81 s    462.08 s   536.9 s     0.2 MB
//
// Two properties are worth stating because the obvious implementation has
// neither:
//
//   - **Resident memory is constant.** A monolithic sieve needs a bit per
//     candidate -- 1.9 GB at P_max = 1.54e10, which is what made the previous
//     `split_primes_up_to` unusable past B ~ 1e10. Segmenting holds one window
//     plus the base primes to sqrt(P_max), which is 0.2 MB at any bound.
//
//   - **The artifact stores (a, b), never p.** Since p = a^2 + b^2 is
//     recoverable, storing it would double the file for nothing: 8 bytes per
//     prime instead of 16, so the full B = 1e12 base is 2.6 GB rather than 5.2.
//     See data/README.md for the on-disk format.
#pragma once

#include "cpu/core/precision.h"
#include "cpu/gaussian/gaussian.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace pte::cpu {

// Segmented sieve of Eratosthenes restricted to split primes.
//
// Yields p = 1 (mod 4), p >= 5, ascending. Does not decompose them -- that is
// FactorBase's job, and separating the two keeps the sieve testable against a
// plain prime list.
class SegmentedSieve {
public:
  // `segment` is the window width in integers. The default keeps the window
  // (one bit per candidate) inside L2 on every machine this targets.
  explicit SegmentedSieve(factor_t p_max, std::size_t segment = 1u << 20);

  factor_t p_max() const { return p_max_; }

  // Calls `cb(factor_t p)` for each split prime <= p_max(), ascending.
  template <typename F> void for_each_split_prime(F &&cb) const;

private:
  factor_t p_max_ = 0;
  std::size_t segment_ = 0;
  std::vector<std::uint32_t> base_; // primes to sqrt(p_max), for the crossing
};

// Split primes with their decompositions, indexed in ascending order of p.
//
// `primes[i]` and `points[i]` describe the same prime, so a DecomposedNorm's
// factor indices address both. Built by sieving, or loaded from an artifact --
// loading skips Cornacchia entirely, which is 86% of the build time.
class FactorBase {
public:
  FactorBase() = default;

  // Sieve and decompose. Throws NotASplitPrime if the sieve ever yields a
  // prime that fails to decompose, which would mean the sieve is wrong.
  static FactorBase Build(factor_t p_max);

  // Read an artifact written by WriteArtifact. The file carries (a, b) pairs
  // only; p is recomputed as a^2 + b^2, so this is a memcpy plus a multiply per
  // prime and no modular arithmetic at all.
  //
  // `p_max` is not stored in the file (it has no header); pass 0 to accept
  // whatever the file contains, or a positive value to truncate to that bound
  // and to assert the file reaches it.
  static FactorBase Load(const std::string &path, factor_t p_max = 0);

  // Sieve, decompose and stream to disk without holding the base in memory.
  // Returns the number of records written.
  static std::size_t WriteArtifact(const std::string &path, factor_t p_max);

  std::size_t size() const { return primes_.size(); }
  bool empty() const { return primes_.empty(); }
  factor_t p_max() const { return p_max_; }

  factor_t prime(factor_index_t i) const { return primes_[i]; }
  const Gauss &point(factor_index_t i) const { return points_[i]; }

  const std::vector<factor_t> &primes() const { return primes_; }
  const std::vector<Gauss> &points() const { return points_; }

private:
  factor_t p_max_ = 0;
  std::vector<factor_t> primes_;
  std::vector<Gauss> points_;
};

// ---------------------------------------------------------------------------
// Implementation of the templated member
// ---------------------------------------------------------------------------

template <typename F> void SegmentedSieve::for_each_split_prime(F &&cb) const {
  if (p_max_ < 5)
    return;

  std::vector<bool> window(segment_);
  for (factor_t lo = 0; lo <= p_max_; lo += segment_) {
    const factor_t hi =
        std::min<factor_t>(lo + static_cast<factor_t>(segment_) - 1, p_max_);
    window.assign(segment_, true);

    for (const std::uint32_t p : base_) {
      const factor_t sq = static_cast<factor_t>(p) * p;
      if (sq > hi)
        break;
      // First multiple of p at or after lo, but never below p^2 -- smaller
      // multiples were crossed out by a smaller base prime.
      factor_t start = ((lo + p - 1) / p) * p;
      if (start < sq)
        start = sq;
      for (factor_t j = start; j <= hi; j += p)
        window[static_cast<std::size_t>(j - lo)] = false;
    }

    const factor_t first = lo < 5 ? 5 : lo;
    for (factor_t n = first; n <= hi; ++n) {
      if (n % 4 == 1 && window[static_cast<std::size_t>(n - lo)])
        cb(n);
    }
  }
}

} // namespace pte::cpu
