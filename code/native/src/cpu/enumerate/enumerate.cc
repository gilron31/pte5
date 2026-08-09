#include "cpu/enumerate/enumerate.h"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace pte::cpu {
namespace {

// GILRO: norm_t and Norm is confusing, we should avoid it by giving better
// names.
// Exact key of a pair. exact_key_t holds 2*L1^2 for any L1 below 9.2e18, which
// is nine million times the design bound; see precision.h.
exact_key_t pair_key(norm_t x, norm_t y) {
  return static_cast<exact_key_t>(x) * x + static_cast<exact_key_t>(y) * y;
}

// The largest prime factor an admissible norm can have. With exponent 1 its
// cofactor needs prod(b+1) >= 4, and the smallest such cofactor is 5*13 = 65;
// every other exponent pattern admits only smaller primes.
//
// This divisor is load-bearing and exactly tight: at B = 10^6 the largest prime
// factor actually occurring is 15377 against a cutoff of 15384, the extremal
// norm being 5 * 13 * 15377 = 999505. An off-by-one here silently deletes
// solutions, so enumerate_test.cc checks it by brute force rather than trusting
// the argument above.
constexpr norm_t kSmallestCofactor = 65;

} // namespace

std::vector<std::uint64_t> split_primes_up_to(std::uint64_t limit) {
  std::vector<std::uint64_t> primes;
  // GILRO: reserve
  if (limit < 5)
    return primes;

  std::vector<bool> composite(limit + 1, false);
  for (std::uint64_t i = 2; i * i <= limit; ++i) {
    if (!composite[i]) {
      for (std::uint64_t j = i * i; j <= limit; j += i)
        composite[j] = true;
    }
  }
  for (std::uint64_t i = 5; i <= limit; ++i) {
    if (!composite[i] && i % 4 == 1)
      primes.push_back(i);
  }
  return primes;
}

std::size_t point_count(const std::vector<std::uint64_t> &primes) {
  // prod(b_i + 1) from the multiset, then halve (rounding up) for conjugation.
  std::size_t divisors = 1;
  std::size_t run = 0;
  for (std::size_t i = 0; i < primes.size(); ++i) {
    ++run;
    if (i + 1 == primes.size() || primes[i + 1] != primes[i]) {
      divisors *= run + 1;
      run = 0;
    }
  }
  return (divisors + 1) / 2;
}

Enumerator::Enumerator(norm_t bound) : bound_(bound) {
  if (bound < 5 * 13 * 17)
    return; // below the smallest admissible norm; factor base stays empty
  // GILRO: throw error maybe?

  primes_ =
      split_primes_up_to(static_cast<std::uint64_t>(bound / kSmallestCofactor));

  // The one place Cornacchia runs: once per prime in the factor base, rather
  // than once per prime factor of every norm.
  prime_points_.reserve(primes_.size());
  for (const std::uint64_t p : primes_)
    prime_points_.push_back(decompose_prime(p));
}

// GILRO: why not iterate by index always, pretty sure this funciton should not
// exist.
std::size_t Enumerator::index_of(std::uint64_t p) const {
  const auto it = std::lower_bound(primes_.begin(), primes_.end(), p);
  if (it == primes_.end() || *it != p) {
    throw std::invalid_argument("Enumerator: prime " + std::to_string(p) +
                                " is not in the factor base for bound " +
                                std::to_string(bound_));
  }
  return static_cast<std::size_t>(it - primes_.begin());
}

std::vector<Gauss> Enumerator::point_set(const Norm &norm) const {
  std::vector<Gauss> prime_points;
  prime_points.reserve(norm.primes.size());
  for (const std::uint64_t p : norm.primes)
    prime_points.push_back(prime_points_[index_of(p)]);
  return point_set_from(prime_points);
}

// GILRO: This should be promoted to something more generic like
// <class Func> transform_points(...)
std::vector<norm_t> nu_values(const std::vector<Gauss> &points) {
  std::vector<norm_t> out;
  out.reserve(points.size());
  for (const Gauss &z : points)
    out.push_back(static_cast<norm_t>(z.nu()));
  return out;
}

// GILRO: What makes a "reference" different that "optimized CPU
// ipmlementation?"
std::vector<Quadruple> meet_reference(const std::vector<norm_t> &nu) {
  std::vector<Quadruple> hits;
  const int n = static_cast<int>(nu.size());
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j)
      for (int k = j + 1; k < n; ++k)
        for (int l = k + 1; l < n; ++l) {
          // The three ways to split {i,j,k,l} into two pairs. nu is injective,
          // so any match here is already a set of four distinct points.
          const exact_key_t kij = pair_key(nu[i], nu[j]);
          const exact_key_t kik = pair_key(nu[i], nu[k]);
          const exact_key_t kil = pair_key(nu[i], nu[l]);
          if (kij == pair_key(nu[k], nu[l]))
            hits.push_back({i, j, k, l});
          if (kik == pair_key(nu[j], nu[l]))
            hits.push_back({i, k, j, l});
          if (kil == pair_key(nu[j], nu[k]))
            hits.push_back({i, l, j, k});
        }
  return hits;
}

bool confirm_candidate(const std::vector<norm_t> &nu, const Quadruple &q) {
  const int n = static_cast<int>(nu.size());
  if (q.a < 0 || q.b < 0 || q.c < 0 || q.d < 0)
    return false;
  if (q.a >= n || q.b >= n || q.c >= n || q.d >= n)
    return false;
  if (q.a == q.b || q.a == q.c || q.a == q.d || q.b == q.c || q.b == q.d ||
      q.c == q.d)
    return false;
  return pair_key(nu[q.a], nu[q.b]) == pair_key(nu[q.c], nu[q.d]);
}

bool is_primitive(const std::vector<Gauss> &points,
                  const std::vector<std::uint64_t> &primes,
                  const Quadruple &q) {
  const Gauss *z[4] = {&points[q.a], &points[q.b], &points[q.c], &points[q.d]};
  for (std::size_t i = 0; i < primes.size(); ++i) {
    // Only repeated primes can divide a point, so skip the first copy of each
    // run and test the prime once.
    if (i + 1 >= primes.size() || primes[i + 1] != primes[i])
      continue;
    if (i > 0 && primes[i] == primes[i - 1])
      continue;

    const auto p = static_cast<coord_t>(primes[i]);
    bool all = true;
    for (int t = 0; t < 4 && all; ++t)
      all = z[t]->divisible_by(p);
    if (all)
      return false; // p divides every point: a scaled copy of L1/p^2
  }
  return true;
}

} // namespace pte::cpu
