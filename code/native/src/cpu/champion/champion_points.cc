#include "cpu/champion/champion_points.h"

#include "cpu/bigint/bigint.h"
#include "cpu/gaussian/gaussian.h"

#include <algorithm>
#include <cstddef>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace pte::cpu {
namespace {

// A Gaussian integer at arbitrary precision. Only the three operations the
// construction needs; see the header for why this is not GaussT<mpz_class>.
struct BigGauss {
  mpz_class re, im;
};

// Both z*q and z*conj(q) from four multiplications rather than eight -- the
// same fused trick as GaussT::mul_and_conj_mul, and it matters more here, where
// a multiplication is a GMP call rather than an instruction.
void MulAndConjMul(const BigGauss &z, const BigGauss &q, BigGauss &out,
                   BigGauss &out_conj) {
  const mpz_class x = z.re * q.re;
  const mpz_class y = z.re * q.im;
  const mpz_class w = z.im * q.re;
  const mpz_class v = z.im * q.im;
  out.re = x - v;
  out.im = y + w;
  out_conj.re = x + v;
  out_conj.im = w - y;
}

// The first-octant representative, m >= n >= 0. This is what makes nu
// injective on the set, which is what licenses the absence of a disjointness
// test in every meet -- see gaussian.h's canonize for the full argument.
BigGauss Canonize(const BigGauss &z) {
  mpz_class a = abs(z.re), b = abs(z.im);
  if (a < b)
    std::swap(a, b);
  return BigGauss{a, b};
}

bool Less(const BigGauss &a, const BigGauss &b) {
  if (a.re != b.re)
    return a.re < b.re;
  return a.im < b.im;
}

bool Equal(const BigGauss &a, const BigGauss &b) {
  return a.re == b.re && a.im == b.im;
}

// D = 2^(k-1) for k distinct primes, and gpu/champion packs a point index into
// 20 bits, so 21 primes is the ceiling that keeps the two consistent.
constexpr std::size_t kMaxPrimes = 21;

} // namespace

std::vector<factor_t> smallest_split_primes(int count) {
  std::vector<factor_t> out;
  if (count <= 0)
    return out;
  for (factor_t p = 5; out.size() < static_cast<std::size_t>(count); p += 4) {
    // Split primes are p = 1 mod 4, so stepping by 4 from 5 visits exactly the
    // right residue class and is_prime does the rest.
    if (is_prime(p))
      out.push_back(p);
  }
  return out;
}

ChampionNorm champion_norm(const std::vector<factor_t> &primes) {
  if (primes.size() > kMaxPrimes)
    throw std::invalid_argument(
        "champion_norm: at most " + std::to_string(kMaxPrimes) +
        " primes; beyond that D exceeds what a point index can address");
  const std::set<factor_t> distinct(primes.begin(), primes.end());
  if (distinct.size() != primes.size())
    throw std::invalid_argument(
        "champion_norm: primes must be distinct; a repeated prime needs the "
        "exponent-sweeping construction in cpu/gaussian, and would break "
        "D = 2^(k-1)");

  ChampionNorm out;
  out.primes.assign(distinct.begin(), distinct.end());
  out.l1 = 1;

  std::vector<BigGauss> points{BigGauss{1, 0}};
  points.reserve(std::size_t{1} << std::min<std::size_t>(primes.size(), 20));
  for (const factor_t p : out.primes) {
    // decompose_prime throws NotASplitPrime on anything that is not a prime
    // congruent to 1 mod 4, which is the loud failure the campaign's audit
    // insisted on: a silent skip here would delete points and the meet would
    // report "no collisions" for a norm it never really built.
    const Gauss q32 = decompose_prime(p);
    const BigGauss q{q32.re, q32.im};
    out.l1 *= p;

    std::vector<BigGauss> next;
    next.reserve(points.size() * 2);
    for (const BigGauss &z : points) {
      BigGauss a, b;
      MulAndConjMul(z, q, a, b);
      next.push_back(std::move(a));
      next.push_back(std::move(b));
    }
    points.swap(next);
  }

  for (BigGauss &z : points)
    z = Canonize(z);
  std::sort(points.begin(), points.end(), Less);
  points.erase(std::unique(points.begin(), points.end(), Equal), points.end());

  out.nu.reserve(points.size());
  for (const BigGauss &z : points) {
    // nu = m^2 - n^2 >= 0 on a canonized point (m >= n >= 0), so the unsigned
    // fixed-width type is the right target; to_uint throws if it does not fit,
    // which is the overflow this whole file exists to make impossible to miss.
    out.nu.push_back(to_uint<4>(z.re * z.re - z.im * z.im));
  }
  return out;
}

} // namespace pte::cpu
