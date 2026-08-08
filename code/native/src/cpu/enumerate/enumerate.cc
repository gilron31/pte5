#include "cpu/enumerate/enumerate.h"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace pte::cpu {
namespace {

// __int128 is a GNU extension; the opt-in keeps -Wpedantic quiet.
__extension__ using i128 = __int128;

// Exact key of a pair, wide enough for any L1 we can enumerate.
i128 pair_key(std::int64_t x, std::int64_t y) {
  return static_cast<i128>(x) * x + static_cast<i128>(y) * y;
}

// True if the rational prime p divides z. A prime dividing *every* point of a
// candidate is what makes it a scaled duplicate; see is_primitive.
bool divides(std::uint64_t p, const Gauss &z) {
  const std::int64_t q = static_cast<std::int64_t>(p);
  return z.re % q == 0 && z.im % q == 0;
}

} // namespace

std::vector<std::uint64_t> split_primes_up_to(std::uint64_t limit) {
  std::vector<std::uint64_t> primes;
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

void for_each_admissible_norm(std::int64_t bound,
                              const std::function<void(const Norm &)> &cb) {
  // The largest prime factor an admissible norm can have. With exponent 1 its
  // cofactor needs prod(b+1) >= 4, and the smallest such cofactor is 5*13 = 65;
  // every other exponent pattern admits only smaller primes.
  if (bound < 5 * 13 * 17)
    return;
  const auto base = split_primes_up_to(static_cast<std::uint64_t>(bound / 65));
  if (base.empty())
    return;

  Norm norm;
  norm.primes.reserve(16);

  // Walks primes in ascending index order, choosing an exponent for each, so
  // every multiset of split primes with product <= bound is visited once.
  const auto dfs = [&](auto &&self, std::size_t start, std::int64_t value,
                       std::size_t divisors) -> void {
    for (std::size_t i = start; i < base.size(); ++i) {
      const std::int64_t p = static_cast<std::int64_t>(base[i]);
      if (value > bound / p)
        break; // primes ascend, so nothing further fits either

      std::int64_t v = value * p;
      std::size_t exponent = 1;
      while (true) {
        norm.primes.push_back(base[i]);
        const std::size_t d = divisors * (exponent + 1);
        if (d >= 8) { // >= 4 points, the minimum for an E_4
          norm.l1 = v;
          cb(norm);
        }
        self(self, i + 1, v, d);

        if (v > bound / p)
          break;
        v *= p;
        ++exponent;
      }
      norm.primes.resize(norm.primes.size() - exponent);
    }
  };
  dfs(dfs, 0, 1, 1);
}

std::vector<std::int64_t> nu_values(const std::vector<Gauss> &points) {
  std::vector<std::int64_t> out;
  out.reserve(points.size());
  for (const Gauss &z : points)
    out.push_back(nu(z));
  return out;
}

std::vector<Quadruple> meet_reference(const std::vector<std::int64_t> &nu) {
  std::vector<Quadruple> hits;
  const int n = static_cast<int>(nu.size());
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j)
      for (int k = j + 1; k < n; ++k)
        for (int l = k + 1; l < n; ++l) {
          // The three ways to split {i,j,k,l} into two pairs. nu is injective,
          // so any match here is already a set of four distinct points.
          const i128 kij = pair_key(nu[i], nu[j]);
          const i128 kik = pair_key(nu[i], nu[k]);
          const i128 kil = pair_key(nu[i], nu[l]);
          if (kij == pair_key(nu[k], nu[l]))
            hits.push_back({i, j, k, l});
          if (kik == pair_key(nu[j], nu[l]))
            hits.push_back({i, k, j, l});
          if (kil == pair_key(nu[j], nu[k]))
            hits.push_back({i, l, j, k});
        }
  return hits;
}

bool verify_exact(const std::vector<std::int64_t> &nu, const Quadruple &q) {
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

    bool all = true;
    for (int t = 0; t < 4 && all; ++t)
      all = divides(primes[i], *z[t]);
    if (all)
      return false; // p divides every point: a scaled copy of L1/p^2
  }
  return true;
}

} // namespace pte::cpu
