#include "cpu/gaussian/gaussian.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace pte::cpu {
namespace {

// GILRO: it is weird to declare this in the .cc file when it should be
// "elevated" to a higher status, it is this type that the precision of the
// entire package relies on. It should be clear.
using u64 = std::uint64_t;
// __int128 is a GNU extension; the opt-in keeps -Wpedantic quiet.
__extension__ using u128 = unsigned __int128;

// GILRO: This is nice, but someday we will want to use larger numbers, guess
// what we will use then (GMP).
u64 mulmod(u64 a, u64 b, u64 m) {
  return static_cast<u64>((static_cast<u128>(a) * b) % m);
}

// Floor of the square root of n. Seeded from the floating-point value, then
// corrected, so it is exact for every 64-bit input.
// GILRO: isn't there a GMP function for this? let's don't reinvent the wheel
// (in a very inefficient manner).
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

// A square root of -1 mod p, for p = 1 (mod 4). Finds a quadratic non-residue n
// by trial and returns n^((p-1)/4), whose square is n^((p-1)/2) = -1 by Euler's
// criterion. Returns 0 if no non-residue is found below the (never reached in
// practice) search bound.
u64 sqrt_minus_one(u64 p) {
  // GILRO: magic number is unacceptable.
  for (u64 n = 2; n < 1000; ++n) {
    if (powmod(n, (p - 1) / 2, p) == p - 1)
      return powmod(n, (p - 1) / 4, p);
  }
  return 0;
}

} // namespace

// GILRO: we are in C++, return std::optional,no, better, throw an error if it
// doesn't decompose and enforce the calling function to check it be itself.
// return a Gauss
bool decompose_prime(u64 p, Gauss *out) {
  if (p % 4 != 1 || p < 5)
    return false;

  const u64 x = sqrt_minus_one(p);
  if (x == 0)
    return false;

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
  if (static_cast<u128>(d) * d != c)
    return false;

  const std::int64_t hi = static_cast<std::int64_t>(std::max(b, d));
  const std::int64_t lo = static_cast<std::int64_t>(std::min(b, d));
  if (lo == 0)
    return false;
  *out = {hi, lo};
  return true;
}

// GILRO: this is weird. this is way to "robust" and memoryless. We compute the
// entire prime decompositions for every norm, very weird. It should be handled
// by a "Enumerator" class or something. SUggest a refactor together with
// enumerate/
std::vector<Gauss> point_set(const std::vector<u64> &primes) {
  std::vector<Gauss> points{{1, 0}};
  points.reserve(std::size_t{1} << std::min<std::size_t>(primes.size(), 20));

  for (const u64 p : primes) {
    Gauss q;
    if (!decompose_prime(p, &q))
      return {}; // caller passed a non-split prime; no points of this norm

    // One pass over the current set, branching each element into z*q and
    // z*conj(q). Iterating over *each copy* of a repeated prime is what sweeps
    // the full exponent range.
    std::vector<Gauss> next;
    next.reserve(points.size() * 2);
    for (const Gauss &z : points) {
      Gauss a;
      Gauss b;
      mul_and_conj_mul(z, q, &a, &b);
      next.push_back(a);
      next.push_back(b);
    }
    points.swap(next);
  }

  for (Gauss &z : points)
    z = canonize(z);
  std::sort(points.begin(), points.end());
  points.erase(std::unique(points.begin(), points.end()), points.end());
  return points;
}

} // namespace pte::cpu
