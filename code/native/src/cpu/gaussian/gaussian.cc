#include "cpu/gaussian/gaussian.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace pte::cpu {
namespace {

using u64 = std::uint64_t;
__extension__ using u128 = unsigned __int128;

u64 mulmod(u64 a, u64 b, u64 m) {
  return static_cast<u64>((static_cast<u128>(a) * b) % m);
}

u64 powmod(u64 base, u64 exp, u64 m) {
  u64 result = 1;
  base %= m;
  while (exp > 0) {
    if (exp & 1)
      result = mulmod(result, base, m);
    base = mulmod(base, base, m);
    exp >>= 1;
  }
  return result;
}

} // namespace

// Seeded from the floating-point square root, then corrected, so it is exact
// for every 64-bit input. Both correction loops are needed: the seed can land
// either side of the true value once n exceeds 2^53 and the double cannot
// represent it exactly.
u64 isqrt(u64 n) {
  if (n == 0)
    return 0;
  u64 r = static_cast<u64>(std::sqrt(static_cast<double>(n)));
  while (r > 0 && static_cast<u128>(r) * r > n)
    --r;
  while (static_cast<u128>(r + 1) * (r + 1) <= n)
    ++r;
  return r;
}

// Deterministic Miller-Rabin.
//
// The witness sets are the standard proven ones: {2,3,5,7} decides every
// n < 3,215,031,751 and {2,3,5,7,11,13} every n < 3,474,749,660,383. The design
// bound puts the factor base at 10^12 / 65 = 1.5e10, comfortably inside the
// second, so this is exact rather than probabilistic for every input the
// enumeration can produce. The twelve-base set covers the rest of the 64-bit
// range, so a caller arriving from outside the enumeration is still safe.
bool is_prime(u64 n) {
  if (n < 2)
    return false;
  for (const u64 small : {2ull, 3ull, 5ull, 7ull, 11ull, 13ull, 17ull, 19ull,
                          23ull, 29ull, 31ull, 37ull}) {
    if (n == small)
      return true;
    if (n % small == 0)
      return false;
  }

  u64 d = n - 1;
  int r = 0;
  while ((d & 1) == 0) {
    d >>= 1;
    ++r;
  }

  static constexpr u64 kSmallSet[] = {2, 3, 5, 7};
  static constexpr u64 kMidSet[] = {2, 3, 5, 7, 11, 13};
  static constexpr u64 kFullSet[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
  const u64 *bases = kFullSet;
  std::size_t count = 12;
  if (n < 3215031751ull) {
    bases = kSmallSet;
    count = 4;
  } else if (n < 3474749660383ull) {
    bases = kMidSet;
    count = 6;
  }

  for (std::size_t i = 0; i < count; ++i) {
    u64 x = powmod(bases[i], d, n);
    if (x == 1 || x == n - 1)
      continue;
    bool composite = true;
    for (int j = 1; j < r; ++j) {
      x = mulmod(x, x, n);
      if (x == n - 1) {
        composite = false;
        break;
      }
    }
    if (composite)
      return false;
  }
  return true;
}

namespace {

// How many candidates to try when looking for a quadratic non-residue mod p.
//
// Never reached for a genuine prime p = 1 (mod 4): exactly half the non-zero
// residues are non-residues, so the chance that the first k trials all fail is
// 2^-k, and 2^-64 is beyond any run this code will ever make.
//
// Note that exhausting it is *not* the primality check -- that would be a
// tempting argument and it is wrong. A composite that is itself a product of
// split primes (85 = 5*17, 125 = 5^3, 221 = 13*17) is a sum of two squares, so
// a non-residue exists, the descent closes, and the result is a well-formed
// decomposition of a non-prime. Those are exactly the composites this code is
// most likely to be handed by mistake, which is why is_prime() above is a real
// Miller-Rabin and not this loop.
constexpr int kMaxNonResidueTrials = 64;

// A square root of -1 mod p, for p = 1 (mod 4). Finds a quadratic non-residue n
// by trial and returns n^((p-1)/4), whose square is n^((p-1)/2) = -1 by Euler's
// criterion. Returns 0 only if p is composite; see above.
u64 sqrt_minus_one(u64 p) {
  for (int t = 0; t < kMaxNonResidueTrials; ++t) {
    const u64 n = static_cast<u64>(t) + 2;
    if (powmod(n, (p - 1) / 2, p) == p - 1)
      return powmod(n, (p - 1) / 4, p);
  }
  return 0;
}

[[noreturn]] void Reject(u64 p, const char *why) {
  throw NotASplitPrime("decompose_prime(" + std::to_string(p) + "): " + why);
}

} // namespace

Gauss decompose_prime(u64 p) {
  if (p < 5 || p % 4 != 1)
    Reject(p, "not congruent to 1 mod 4 (2 is ramified, 3 mod 4 is inert)");
  if (!is_prime(p))
    Reject(p, "composite; a product of split primes also has a decomposition, "
              "so the descent below would have returned one");

  const u64 x = sqrt_minus_one(p);
  if (x == 0)
    Reject(p, "no quadratic non-residue found in 64 trials, so p is composite");

  // Cornacchia: run the Euclidean algorithm on (p, x) and stop at the first
  // remainder below sqrt(p). Because x^2 = -1 (mod p), the pair (a, b) of the
  // last two remainders then satisfies a^2 + b^2 = p.
  u64 a = p;
  u64 b = x;
  while (static_cast<u128>(b) * b > p) {
    const u64 r = a % b;
    a = b;
    b = r;
  }

  const u64 c = static_cast<u64>(p - static_cast<u128>(b) * b);
  const u64 d = isqrt(c);
  if (static_cast<u128>(d) * d != c || d == 0)
    Reject(p, "Cornacchia's descent did not close, so p is composite");

  const auto hi = static_cast<coord_t>(std::max(b, d));
  const auto lo = static_cast<coord_t>(std::min(b, d));
  return {hi, lo};
}

std::vector<Gauss> point_set_from(const std::vector<Gauss> &prime_points) {
  std::vector<Gauss> points{{1, 0}};
  points.reserve(std::size_t{1}
                 << std::min<std::size_t>(prime_points.size(), 20));

  for (const Gauss &q : prime_points) {
    // One pass over the current set, branching each element into z*q and
    // z*conj(q). Iterating over *each copy* of a repeated prime is what sweeps
    // the full exponent range.
    std::vector<Gauss> next;
    next.reserve(points.size() * 2);
    for (const Gauss &z : points) {
      Gauss a;
      Gauss b;
      z.mul_and_conj_mul(q, &a, &b);
      next.push_back(a);
      next.push_back(b);
    }
    points.swap(next);
  }

  for (Gauss &z : points)
    z = z.canonize();
  std::sort(points.begin(), points.end());
  points.erase(std::unique(points.begin(), points.end()), points.end());
  return points;
}

std::vector<Gauss> point_set(const std::vector<u64> &primes) {
  std::vector<Gauss> prime_points;
  prime_points.reserve(primes.size());
  for (const u64 p : primes)
    prime_points.push_back(decompose_prime(p));
  return point_set_from(prime_points);
}

} // namespace pte::cpu
