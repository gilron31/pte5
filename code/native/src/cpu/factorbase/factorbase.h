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
// Three properties worth stating because the obvious implementation has none of
// them:
//
//   - **Sieving is segmented**, so resident memory during construction is
//     constant. A monolithic sieve needs a bit per candidate -- 1.9 GB at
//     P_max = 1.54e10, which is what made the previous `split_primes_up_to`
//     unusable past B ~ 1e10.
//
//   - **The prime is never stored.** Only (a, b) is; p = a^2 + b^2 is recovered
//     with two multiplies. That is 8 bytes per prime rather than 24, and it is
//     the difference between 8.2 GB and 2.7 GB at the design bound. It also
//     *improves* locality: prime and point now share a cache line instead of
//     living in two arrays a gigabyte apart.
//
//   - **The artifact can be mapped rather than read.** Map() mmaps the file, so
//     the base costs no heap at all -- pages are demand-loaded, evictable under
//     pressure, and shared between concurrent runs. At B = 10^12 this is the
//     difference between a 2.7 GB allocation and none.
//
// Measured resident set for the whole B = 10^10 campaign: 141 MB with Load,
// where the base is ~104 MB of it.
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
// One 8-byte record per prime, holding (a, b) with p = a^2 + b^2, a >= b > 0.
// The records are either owned (Build, Load) or mapped from an artifact (Map);
// every accessor works the same either way.
class FactorBase {
public:
  // On-disk and in-memory record. See data/README.md for the file format --
  // they are the same bytes, which is what lets Map() hand out pointers
  // straight into the mapping.
  struct Record {
    std::uint32_t a, b;
  };
  static_assert(sizeof(Record) == 8, "artifact record must be exactly 8 bytes");

  FactorBase() = default;
  ~FactorBase();
  FactorBase(FactorBase &&) noexcept;
  FactorBase &operator=(FactorBase &&) noexcept;
  FactorBase(const FactorBase &) = delete;
  FactorBase &operator=(const FactorBase &) = delete;

  // Sieve and decompose. Throws NotASplitPrime if the sieve ever yields a
  // prime that fails to decompose, which would mean the sieve is wrong.
  static FactorBase Build(factor_t p_max);

  // Read an artifact into memory. `p_max` of 0 accepts whatever the file
  // contains; a positive value truncates to that bound and asserts the file
  // reaches it.
  static FactorBase Load(const std::string &path, factor_t p_max = 0);

  // Map an artifact instead of reading it. Same validation as Load on the
  // prefix actually used, but no heap: at B = 10^12 the base is 2.7 GB of
  // page cache rather than 2.7 GB of anonymous memory.
  //
  // The mapping is read-only and private. The file must outlive the object.
  static FactorBase Map(const std::string &path, factor_t p_max = 0);

  // Sieve, decompose and stream to disk without holding the base in memory.
  // Returns the number of records written.
  static std::size_t WriteArtifact(const std::string &path, factor_t p_max);

  std::size_t size() const { return count_; }
  bool empty() const { return count_ == 0; }
  factor_t p_max() const { return p_max_; }
  bool is_mapped() const { return mapping_ != nullptr; }

  // p = a^2 + b^2, recomputed rather than stored. Two multiplies against a
  // load from a second array that would not be in cache anyway.
  factor_t prime(factor_index_t i) const {
    const Record &r = records_[i];
    return static_cast<factor_t>(r.a) * r.a + static_cast<factor_t>(r.b) * r.b;
  }
  Gauss point(factor_index_t i) const {
    const Record &r = records_[i];
    return Gauss{static_cast<coord_t>(r.a), static_cast<coord_t>(r.b)};
  }

  const Record *records() const { return records_; }

  // Exhaustive structural check: every record in order, a >= b > 0. Faults in
  // the entire mapping when mapped, so it is a deliberate audit step rather
  // than something Map does for you.
  void validate_all(const std::string &path = "<factor base>") const;

  // Materializes every prime into a vector. O(size()) memory -- 2.7 GB at the
  // design bound -- so this is for tests and small bounds, not for a campaign.
  std::vector<factor_t> primes_vector() const;

  // Index of the first prime greater than `limit`, i.e. how many primes are at
  // or below it. Binary search over the recomputed primes.
  std::size_t upper_bound(factor_t limit) const;

private:
  void Reset();
  void ValidateSample(const std::string &path) const;

  const Record *records_ = nullptr; // owned_.data(), or into the mapping
  std::size_t count_ = 0;
  factor_t p_max_ = 0;

  std::vector<Record> owned_;      // Build / Load
  void *mapping_ = nullptr;        // Map
  std::size_t mapping_bytes_ = 0;
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
