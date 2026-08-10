"""Checks for the load-bearing new claims in drafts/gpu_implementation_plan.md.

Run:  python3 drafts/verify_gpu_plan.py

Four claims, in the order the plan makes them:

  1. (plan 2.1)  The S-D tower recursion  T_d = { t^2 - L_d : t in A_{d-1} }
     equals the direct definition "E_d(x) - c splits into linear factors over
     F_p", for every prefix (L_1..L_d).  Checked against explicit polynomial
     root-counting with multiplicity.

  2. (plan 2.1)  The state collapse: distinct reachable T-sets per level, which
     is what turns a p^5 count into an O(p^4/64) one.

  3. (plan 2.2)  THE LANDMINE.  The distinct-root variant of the recursion gives
     zero depth-5 towers for almost every p <= 113 -- which would be an
     unconditional nonexistence proof if the model were sound.  It is not: the
     16 roots of the section-7 exemplar E_4, distinct over Z, collide mod p at
     more than half of those primes.  An integer solution therefore need not
     reduce to a distinct-root tower, and sigma_p^distinct = 0 excludes nothing.

  4. (plan 3.1)  Windowing the champion meet by key *value* with a two-pointer
     over sorted squares enumerates every pair exactly once, so the cost is
     C(D,2) rather than R * C(D,2).
"""

from __future__ import annotations

import bisect
import itertools
import math
import random
import sys

# The section-7 exemplar E_4, from paper.md / drafts/verify_claims.py.
EXEMPLAR = (
    33065996717,
    489628056848329146064,
    175480010455650701584492675662518592000000,
    40042900368028062136207226327668992000000 ** 2,
)

PRIMES_113 = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
              67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113)


# ---------------------------------------------------------------------------
# 1. The recursion
# ---------------------------------------------------------------------------

def tower_poly(Ls: tuple[int, ...], p: int) -> list[int]:
    """E_d as coefficients over F_p, low order first.  E_0(x) = x."""
    F = [0, 1]
    for L in Ls:
        sq = [0] * (2 * len(F) - 1)
        for i, ai in enumerate(F):
            if ai:
                for j, aj in enumerate(F):
                    sq[i + j] = (sq[i + j] + ai * aj) % p
        F = sq
        F[0] = (F[0] - L) % p
    return F


def splits_with_multiplicity(coeffs: list[int], p: int) -> bool:
    """Does the polynomial factor into linear factors over F_p?"""
    deg = len(coeffs) - 1
    total, c = 0, coeffs[:]
    for r in range(p):
        while len(c) > 1:
            n = len(c) - 1
            q, rem = [0] * n, c[n]
            for i in range(n - 1, -1, -1):
                q[i] = rem
                rem = (c[i] + rem * r) % p
            if rem != 0:
                break
            c, total = q, total + 1
    return total == deg


def brute_T(Ls: tuple[int, ...], p: int, distinct: bool) -> frozenset[int]:
    F = tower_poly(Ls, p)
    deg = len(F) - 1
    out = set()
    for c in range(p):
        if distinct:
            roots = {x for x in range(p)
                     if sum(F[i] * pow(x, i, p) for i in range(len(F))) % p == c}
            if len(roots) == deg:
                out.add(c)
        else:
            g = F[:]
            g[0] = (g[0] - c) % p
            if splits_with_multiplicity(g, p):
                out.add(c)
    return frozenset(out)


def dp_T(Ls: tuple[int, ...], p: int, distinct: bool) -> frozenset[int]:
    """The recursion under test."""
    T = frozenset(range(p))
    for L in Ls:
        A = frozenset(t for t in T
                      if (-t) % p in T and (not distinct or t != 0))
        T = frozenset((t * t - L) % p for t in A)
    return T


def check_recursion() -> bool:
    ok = True
    for distinct in (False, True):
        label = "distinct roots" if distinct else "with multiplicity"
        for p in (3, 5, 7, 11, 13):
            for d in (1, 2, 3):
                if p ** d > 3000:
                    continue
                for Ls in itertools.product(range(p), repeat=d):
                    if dp_T(Ls, p, distinct) != brute_T(Ls, p, distinct):
                        print(f"  MISMATCH p={p} Ls={Ls} ({label})")
                        ok = False
        print(f"  recursion == brute force, {label:18s} "
              f"p = 3..13, d = 1,2,3")
    return ok


# ---------------------------------------------------------------------------
# 2 + 3. The level DP, and the two local models
# ---------------------------------------------------------------------------

def level_dp(p: int, depth: int, distinct: bool):
    """Returns (tuples_at_each_level, distinct_states_at_each_level)."""
    layer = {frozenset(range(p)): 1}
    counts, states = [], []
    for _ in range(depth):
        states.append(len(layer))
        nxt: dict[frozenset[int], int] = {}
        total = 0
        for T, mult in layer.items():
            A = frozenset(t for t in T
                          if (-t) % p in T and (not distinct or t != 0))
            total += mult * len({(t * t) % p for t in A})
            for L in range(p):
                nT = frozenset((t * t - L) % p for t in A)
                nxt[nT] = nxt.get(nT, 0) + mult
        counts.append(total)
        layer = nxt
    return counts, states


def check_local_models() -> bool:
    print("\n  depth-5 tower counts over F_p, both local models:\n")
    print(f"  {'p':>4}  {'with multiplicity':>18}  {'distinct roots':>15}"
          f"   states/level (distinct)")
    zeros = 0
    for p in PRIMES_113:
        mult, _ = level_dp(p, 5, distinct=False)
        dist, st = level_dp(p, 5, distinct=True)
        if dist[-1] == 0:
            zeros += 1
        print(f"  {p:>4}  {mult[-1]:>18}  {dist[-1]:>15}   {st}")
    print(f"\n  distinct-root model gives ZERO depth-5 towers at "
          f"{zeros}/{len(PRIMES_113)} primes <= 113.")
    print("  Read naively that is an unconditional nonexistence proof. It is not:")
    return check_reduction_collides()


def tower_roots(Ls: tuple[int, ...]) -> list[int]:
    """The 2^d integer roots of a tower, by peeling from the top."""
    vals = [0]
    for L in reversed(Ls):
        nxt = []
        for v in vals:
            t = L + v
            r = math.isqrt(t)
            assert r * r == t, "tower does not peel: not a perfect square"
            nxt += [r, -r] if r else [0]
        vals = nxt
    return vals


def check_reduction_collides() -> bool:
    roots = tower_roots(EXEMPLAR)
    assert len(set(roots)) == 16, "exemplar should have 16 distinct roots over Z"
    bad = [(p, len({r % p for r in roots})) for p in PRIMES_113
           if len({r % p for r in roots}) < 16]
    print(f"  the exemplar E_4 (L1 = {EXEMPLAR[0]}) has 16 roots, distinct over Z,")
    print(f"  but only {', '.join(f'{d} mod {p}' for p, d in bad[:6])}, ...")
    print(f"  -> collides at {len(bad)}/{len(PRIMES_113)} primes <= 113.")
    print("  A real solution fails the distinct-root local test at more than half")
    print("  of them, so sigma_p^distinct = 0 excludes nothing. The nonexistence")
    print("  branch must use the multiplicity-allowing model.")
    return len(bad) > 0


# ---------------------------------------------------------------------------
# 4. Value-windowed two-pointer for the champion meet
# ---------------------------------------------------------------------------

def check_windowing(D: int = 400, R: int = 7, seed: int = 1) -> bool:
    random.seed(seed)
    sq = sorted(random.randrange(10 ** 6) for _ in range(D))
    lo, hi = 2 * sq[0], 2 * sq[-1] + 1
    edges = [lo + (hi - lo) * r // R for r in range(R + 1)]
    seen, work = set(), 0
    for r in range(R):
        a, b = edges[r], edges[r + 1]
        for i in range(D):
            j0 = max(i + 1, bisect.bisect_left(sq, a - sq[i]))
            j1 = bisect.bisect_left(sq, b - sq[i])
            for j in range(j0, j1):
                work += 1
                if (i, j) in seen:      # a pair enumerated twice is wasted work
                    print(f"  DUPLICATE pair {(i, j)}")
                    return False
                seen.add((i, j))
    allpairs = {(i, j) for i in range(D) for j in range(i + 1, D)}
    if seen != allpairs:               # a pair missed is a LOST SOLUTION
        print(f"  MISSED {len(allpairs - seen)} pairs")
        return False
    print(f"  D = {D}, R = {R}: enumerated {work} pairs, all distinct, "
          f"covering all C(D,2) = {len(allpairs)}")
    print(f"  hash-windowing would evaluate {R * len(allpairs)} keys "
          f"-- {R}x more, and R is 10^3 in the sized design")
    return True


def main() -> int:
    print("1. S-D tower recursion (plan 2.1)")
    ok = check_recursion()
    print("\n2/3. Local models and the distinct-root landmine (plan 2.2)")
    ok &= check_local_models()
    print("\n4. Value-windowed champion meet (plan 3.1)")
    ok &= check_windowing()
    print("\nALL CHECKS PASSED" if ok else "\nFAILURES ABOVE")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
