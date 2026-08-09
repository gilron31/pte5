// Layer 3: recovering an E_4 tower from a meet hit, and proving it.
//
// Everything upstream of this file works with nu values and 64-bit keys, and a
// hit there is a *candidate*. This is where a candidate becomes a result, with
// a machine-checked certificate attached. paper.md's Appendix carries an
// explicit TODO to confirm the degree-16 identity and the ladder divisibilities
// "before citing", so this is a paper obligation rather than hygiene.
//
// ---------------------------------------------------------------------------
// The peel
// ---------------------------------------------------------------------------
//
// From four canonized points z_k = m_k + n_k i of norm L1, paired as
// {z_0, z_1} | {z_2, z_3}:
//
//   roots      r_{2k}, r_{2k+1} = m_k -+ n_k        (8 of them; r^2 + r'^2 = 2L1)
//   derived    u_k = 2 m_k n_k
//   L2         (u_0^2 + u_1^2) / 2                  ( = (u_2^2 + u_3^2) / 2 )
//   v_j        (u_{2j}^2 - u_{2j+1}^2) / 2
//   L3         (v_0^2 + v_1^2) / 2
//   sqrt(L4)   |(v_0^2 - v_1^2) / 2|
//
// **The meet condition is the L2 condition.** nu(z)^2 + u^2 = (m^2-n^2)^2 +
// (2mn)^2 = (m^2+n^2)^2 = L1^2, so for four points of the same norm
//
//     u_a^2 + u_b^2 = u_c^2 + u_d^2   <=>   nu_a^2 + nu_b^2 = nu_c^2 + nu_d^2,
//
// which is exactly what the meet searches for. That equivalence is why L2 is
// well defined for every hit, and why every hit is genuinely an E_4: L3 and L4
// are then *determined*, with no further condition to satisfy. Nothing here can
// discover that a hit was not a solution -- what it can discover is that this
// code computed the wrong tower, which is the whole point of checking.
//
// ---------------------------------------------------------------------------
// Widths
// ---------------------------------------------------------------------------
//
// u <= L1 <= 10^12 still fits a 64-bit word, but L2 ~ u^2 reaches 10^24 and
// L3 ~ v^2 reaches 10^48, past even __int128. The degree-16 identity's
// coefficients run to ~300 bits. This is the layer cpu/core/precision.h
// reserves GMP for, and the one place in the project where a heap allocation
// per operation is affordable: it runs once per solution, and there are a
// handful.
#pragma once

#include "cpu/core/precision.h"
#include "cpu/core/solution.h"
#include "cpu/gaussian/gaussian.h"

#include <gmpxx.h>

#include <array>
#include <span>
#include <string>

namespace pte::cpu {

// Required divisors of the PTE constant, from [CMSV, Table 1] via paper.md 4.
//
//   C'_8 contains 2^6 * 3^3 * 5^2 * 7^2 * 11 * 13 = 302,702,400
//   C'_4 = 2^2 * 3^2 = 36
//
// Theorem 4.1 divides each by two: C'_{2^k}/2 divides the level-k radical.
inline constexpr long kC8Half = 151'351'200;
inline constexpr long kC4Half = 18;

// A recovered E_4 tower.
struct Tower {
  norm_t l1 = 0;
  std::array<Gauss, 4> points{}; // in {0,1} | {2,3} pairing order
  std::array<mpz_class, 8> roots{}; // m_k -+ n_k, pair k at 2k and 2k+1
  mpz_class l2, l3, sqrt_l4;
};

// The outcome of Theorem 4.1 on one tower, spelled out rather than collapsed to
// a bool, because which clause fails says what went wrong: a bad pairing breaks
// the radicals, a bad sign breaks only the divisibility.
struct Ladder {
  bool sqrt_l4_divisible = false; // C'_8 / 2 | sqrt(L4)
  bool plus_is_square = false;    // L3 + sqrt(L4) is a perfect square...
  bool minus_is_square = false;   // ...and so is L3 - sqrt(L4)
  bool plus_divisible = false;    // C'_4 / 2 | sqrt(L3 + sqrt(L4))
  bool minus_divisible = false;

  bool ok() const {
    return sqrt_l4_divisible && plus_is_square && minus_is_square &&
           plus_divisible && minus_divisible;
  }
};

// Everything needed to publish a solution.
struct Certificate {
  Tower tower;
  bool identity_holds = false; // the degree-16 identity, coefficient by
                               // coefficient
  Ladder ladder;

  bool valid() const { return identity_holds && ladder.ok(); }

  // Human-readable account of what failed, empty when valid.
  std::string failure() const;
};

// Recover the tower from four points read as {0,1} | {2,3}.
//
// Throws std::invalid_argument if the four norms disagree, if a point is not
// canonized, or if the meet condition fails. That last check is the exact
// arithmetic confirmation of a filter hit -- it subsumes confirm_candidate, in
// GMP rather than __int128 -- so a truncated-key false positive cannot reach a
// certificate.
Tower recover_tower(const std::array<Gauss, 4> &points);

// Expand prod(x^2 - r_k^2) and (((x^2 - L1)^2 - L2)^2 - L3)^2 - L4 as degree-16
// polynomials over Z and compare all seventeen coefficients.
bool verify_identity(const Tower &t);

// Theorem 4.1 on this tower.
Ladder verify_ladder(const Tower &t);

// Recover and prove, in one call.
Certificate certify(const std::array<Gauss, 4> &points);

// Pull a solution's four points out of a point set, in pairing order. Bounds
// checked, so a garbled hit from a buggy kernel throws rather than reading
// past the end.
std::array<Gauss, 4> quadruple(std::span<const Gauss> points, const E4Sol &q);

} // namespace pte::cpu
