"""Probe for G1 -- the prototype behind drafts/local_model_plan.md section 2.

Run:  python3 drafts/probe_local_model.py

Does the tower DP mod p^k force node valuations to grow with k?

The plan's DP is mod p.  This runs the identical recursion mod p^k and adds one
thing: a cap on the p-adic valuation of every node value (every square root
taken).  A genuine integer tower has finite node valuations, so if the minimum
achievable cap grows without bound in k, every p-adic tower is degenerate and no
integer tower exists.  If it stabilises, that stabilised value is a forced
divisor -- the divisibility ladder, computed instead of proved.
"""

from collections import defaultdict
import sys


def dp(p, k, depth, vcap=None):
    """Returns per-level (#L-tuples with a split tower, #distinct states).

    vcap: if not None, every node value t must satisfy v_p(t) <= vcap.
    """
    m = p**k

    def val(t):
        if t % m == 0:
            return k  # indistinguishable from infinity at this precision
        v = 0
        while t % p == 0:
            t //= p
            v += 1
        return v

    allowed = [t for t in range(m) if vcap is None or val(t) <= vcap]
    # squares of allowed values, as a lookup from t to t^2
    sq = {t: (t * t) % m for t in allowed}

    full = frozenset(range(m))
    layer = {full: 1}
    counts, states = [], []
    for _ in range(depth):
        states.append(len(layer))
        nxt = defaultdict(int)
        total = 0
        for T, mult in layer.items():
            A = [t for t in allowed if t in T and (-t) % m in T]
            images = {sq[t] for t in A}
            total += mult * len(images)  # L_d completing this prefix
            for L in range(m):
                nT = frozenset((s - L) % m for s in images)
                nxt[nT] += mult
        counts.append(total)
        layer = dict(nxt)
    return counts, states


def min_vcap(p, k, depth):
    """Smallest cap admitting a tower mod p^k; k means 'none below k'."""
    for v in range(k + 1):
        c, _ = dp(p, k, depth, vcap=v)
        if c[-1] > 0:
            return v, c[-1]
    return None, 0


# The section-7 exemplar E_4: a *real* depth-4 tower, so it must appear in the
# table above once k exceeds its node valuations.  Printing them is the control
# that makes a null result at depth 5 mean anything -- see the plan's section 2.
EXEMPLAR = (
    33065996717,
    489628056848329146064,
    175480010455650701584492675662518592000000,
    40042900368028062136207226327668992000000**2,
)


def exemplar_valuations(primes=(2, 3, 5, 7)):
    import math

    levels, vals = [], [0]
    for L in reversed(EXEMPLAR):
        nxt = []
        for x in vals:
            t = L + x
            r = math.isqrt(t)
            assert r * r == t, "tower does not peel: not a perfect square"
            nxt += [r, -r] if r else [0]
        vals = nxt
        levels.append(list(nxt))

    def v(n, p):
        if n == 0:
            return None
        c = 0
        while n % p == 0:
            n //= p
            c += 1
        return c

    print("exemplar E_4 node valuations (level 1 = roots, level 4 = sqrt(L_4)):")
    for p in primes:
        prof = [sorted({v(abs(x), p) for x in lv}) for lv in levels]
        top = max(max(s) for s in prof)
        print(
            f"  p={p:<3} "
            + "  ".join(f"L{4-i}:{s}" for i, s in enumerate(prof))
            + f"   -> ht_{p} = {top}, so precision k > {top} is required"
        )


def main():
    exemplar_valuations()
    print()
    print("p^k  depth   towers(uncapped)   min node-valuation cap admitting one")
    for p in (3, 5, 7):
        for k in (1, 2, 3):
            if p**k > 200:
                continue
            for depth in (3, 4, 5):
                counts, states = dp(p, k, depth)
                v, n = min_vcap(p, k, depth)
                print(
                    f"{p}^{k}={p**k:<4} d={depth}   {counts[-1]:>14}   "
                    f"vcap={v} (towers {n}), states/level {states}"
                )
            print()


if __name__ == "__main__":
    sys.exit(main())
