#include "cpu/confirm/confirm.h"

#include <stdexcept>
#include <string>
#include <vector>

namespace pte::cpu {
namespace {

// A dense polynomial over Z, coefficient i multiplying x^i. Degree 16 with
// ~300-bit coefficients, built once per solution -- so this is written for
// clarity, not for speed.
using Poly = std::vector<mpz_class>;

Poly Mul(const Poly &a, const Poly &b) {
  Poly out(a.size() + b.size() - 1);
  for (std::size_t i = 0; i < a.size(); ++i)
    for (std::size_t j = 0; j < b.size(); ++j)
      out[i + j] += a[i] * b[j];
  return out;
}

// p -> p^2 - c, the tower's one structural operation.
Poly SquareMinus(const Poly &p, const mpz_class &c) {
  Poly out = Mul(p, p);
  out[0] -= c;
  return out;
}

// A perfect-square test that also hands back the root. mpz_sqrt truncates, so
// squaring back is what makes this exact rather than approximate.
bool ExactSqrt(const mpz_class &n, mpz_class &root) {
  if (sgn(n) < 0)
    return false;
  mpz_sqrt(root.get_mpz_t(), n.get_mpz_t());
  return root * root == n;
}

} // namespace

std::array<Gauss, 4> quadruple(std::span<const Gauss> points, const E4Sol &q) {
  const int idx[4] = {q.a, q.b, q.c, q.d};
  std::array<Gauss, 4> out{};
  for (int k = 0; k < 4; ++k) {
    if (idx[k] < 0 || static_cast<std::size_t>(idx[k]) >= points.size()) {
      throw std::invalid_argument("quadruple: index " +
                                  std::to_string(idx[k]) +
                                  " is outside a point set of size " +
                                  std::to_string(points.size()));
    }
    out[static_cast<std::size_t>(k)] = points[static_cast<std::size_t>(idx[k])];
  }
  return out;
}

Tower recover_tower(const std::array<Gauss, 4> &points) {
  Tower t;
  t.points = points;

  // All four must share a norm, and be canonized: the roots are m - n and
  // m + n, and a point with n > m would hand back a negative first root. The
  // identity would survive it (only r^2 appears) but the tower would no longer
  // be the one paper.md 5 describes.
  const auto l1 = static_cast<norm_t>(points[0].norm());
  for (int k = 0; k < 4; ++k) {
    const Gauss &z = points[static_cast<std::size_t>(k)];
    if (static_cast<norm_t>(z.norm()) != l1) {
      throw std::invalid_argument(
          "recover_tower: point " + std::to_string(k) + " has norm " +
          std::to_string(static_cast<long long>(z.norm())) + ", not " +
          std::to_string(l1));
    }
    if (z.re < z.im || z.im < 0) {
      throw std::invalid_argument("recover_tower: point " + std::to_string(k) +
                                  " is not canonized (m >= n >= 0)");
    }
  }
  t.l1 = l1;

  // u_k = 2 m_k n_k, and the eight roots m_k -+ n_k.
  std::array<mpz_class, 4> u;
  for (int k = 0; k < 4; ++k) {
    const Gauss &z = points[static_cast<std::size_t>(k)];
    u[static_cast<std::size_t>(k)] = mpz_class(2) * z.re * z.im;
    t.roots[static_cast<std::size_t>(2 * k)] = mpz_class(z.re) - z.im;
    t.roots[static_cast<std::size_t>(2 * k + 1)] = mpz_class(z.re) + z.im;
  }

  // The meet condition, exactly. Above kExactKeyBound the enumeration compares
  // keys modulo 2^64, so this is where a wraparound false positive dies.
  const mpz_class left = u[0] * u[0] + u[1] * u[1];
  const mpz_class right = u[2] * u[2] + u[3] * u[3];
  if (left != right) {
    throw std::invalid_argument(
        "recover_tower: the pairs do not meet (u_0^2 + u_1^2 != u_2^2 + u_3^2 "
        "at L1 = " +
        std::to_string(l1) + "), so there is no L2");
  }

  // Every division below is exact: u is even, so u^2 = 0 mod 4 and each
  // half-difference is an integer.
  t.l2 = left / 2;
  const mpz_class v0 = (u[0] * u[0] - u[1] * u[1]) / 2;
  const mpz_class v1 = (u[2] * u[2] - u[3] * u[3]) / 2;
  t.l3 = (v0 * v0 + v1 * v1) / 2;
  t.sqrt_l4 = abs((v0 * v0 - v1 * v1) / 2);
  return t;
}

bool verify_identity(const Tower &t) {
  // prod (x^2 - r_k^2)
  Poly lhs{mpz_class(1)};
  for (const mpz_class &r : t.roots) {
    const Poly factor{-(r * r), mpz_class(0), mpz_class(1)};
    lhs = Mul(lhs, factor);
  }

  // (((x^2 - L1)^2 - L2)^2 - L3)^2 - L4
  Poly rhs{-mpz_class(t.l1), mpz_class(0), mpz_class(1)};
  rhs = SquareMinus(rhs, t.l2);
  rhs = SquareMinus(rhs, t.l3);
  rhs = SquareMinus(rhs, t.sqrt_l4 * t.sqrt_l4);

  if (lhs.size() != rhs.size())
    return false;
  for (std::size_t i = 0; i < lhs.size(); ++i)
    if (lhs[i] != rhs[i])
      return false;
  return true;
}

Ladder verify_ladder(const Tower &t) {
  Ladder out;
  out.sqrt_l4_divisible = mpz_divisible_ui_p(t.sqrt_l4.get_mpz_t(),
                                             static_cast<unsigned long>(kC8Half)) != 0;

  // The level-2 radicals. Both must be perfect squares -- the tree requires it,
  // and their roots are |v_0| and |v_1| -- so taking the square root here is a
  // real check on l3 and sqrt_l4 rather than a restatement of how they were
  // built.
  mpz_class root;
  const mpz_class plus = t.l3 + t.sqrt_l4;
  if (ExactSqrt(plus, root)) {
    out.plus_is_square = true;
    out.plus_divisible =
        mpz_divisible_ui_p(root.get_mpz_t(),
                           static_cast<unsigned long>(kC4Half)) != 0;
  }
  const mpz_class minus = t.l3 - t.sqrt_l4;
  if (ExactSqrt(minus, root)) {
    out.minus_is_square = true;
    out.minus_divisible =
        mpz_divisible_ui_p(root.get_mpz_t(),
                           static_cast<unsigned long>(kC4Half)) != 0;
  }
  return out;
}

Certificate certify(const std::array<Gauss, 4> &points) {
  Certificate c;
  c.tower = recover_tower(points);
  c.identity_holds = verify_identity(c.tower);
  c.ladder = verify_ladder(c.tower);
  return c;
}

std::string Certificate::failure() const {
  if (valid())
    return {};
  std::string why;
  const auto note = [&why](const char *s) {
    if (!why.empty())
      why += "; ";
    why += s;
  };
  if (!identity_holds)
    note("the degree-16 identity does not hold");
  if (!ladder.sqrt_l4_divisible)
    note("C'_8/2 does not divide sqrt(L4)");
  if (!ladder.plus_is_square)
    note("L3 + sqrt(L4) is not a perfect square");
  if (!ladder.minus_is_square)
    note("L3 - sqrt(L4) is not a perfect square");
  if (ladder.plus_is_square && !ladder.plus_divisible)
    note("C'_4/2 does not divide sqrt(L3 + sqrt(L4))");
  if (ladder.minus_is_square && !ladder.minus_divisible)
    note("C'_4/2 does not divide sqrt(L3 - sqrt(L4))");
  return "L1 = " + std::to_string(tower.l1) + ": " + why;
}

} // namespace pte::cpu
