"""Upper bounds on the p-adic valuation floor at depth 5, from explicit Z_p-towers.

Construction: an E_4 seed with tree-ordered roots r' (either catalog form),
a scale mu = p^v and a bottom value L such that every L +- mu*r'_k is a
square in Z_p.  Then the E_5 with leaves sqrt(L +- mu r'_k) is a nondegenerate
Z_p-tower (Corollary 3.2 over Z_p), with node valuations
    level 1: v(L +- mu r')/2,   level j >= 2: 2^(j-2) v(mu) + v(seed level j-1).
L = 1, v(mu) = 1 always works for p odd (and v(mu) = 3 for p = 2); the search
below tries v(mu) = 0..6 and L = 1..63 and keeps the smallest max valuation.
Valuations are exact (Fractions): the point form has half-integer nodes.
"""
import csv
from fractions import Fraction
import sys
from math import isqrt


def vp(n, p):
    n = Fraction(n)
    if n.denominator != 1:
        return vp(n.numerator, p) - vp(n.denominator, p)
    n = n.numerator
    if n == 0:
        return 10 ** 9
    n, v = abs(n), 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def derived(roots):
    return [abs(Fraction(roots[2 * k]) ** 2 - Fraction(roots[2 * k + 1]) ** 2) / 2
            for k in range(len(roots) // 2)]


def seed_levels(roots):
    levels = [roots]
    while len(levels[-1]) > 1:
        levels.append(derived(levels[-1]))
    return levels  # [roots(8), level2(4), level3(2), level4(1)]


def seeds_from_catalog(path):
    out = []
    for row in csv.DictReader(open(path)):
        L1 = int(row['l1'])
        nus = [int(row[f'nu{i}']) for i in range(4)]
        pts = []
        for nu in nus:
            m, n = isqrt((L1 + nu) // 2), isqrt((L1 - nu) // 2)
            assert m * m + n * n == L1 and m * m - n * n == nu
            pts.append((m, n))
        twisted = [v for m, n in pts for v in (m - n, m + n)]   # roots with r^2+s^2 = 2 L1
        point = [v for m, n in pts for v in (m, n)]             # the 45-degree twin
        out.append((L1, 'twisted', twisted))
        out.append((L1, 'point', point))
    return out


def is_square_Zp(a, p):
    """a in Q (nonzero): is it a square in Q_p?"""
    v = vp(a, p)
    if v % 2:
        return False
    a = Fraction(a)
    u = a / Fraction(p) ** v
    u = u.numerator * pow(u.denominator, -1, 8 if p == 2 else p)
    if p == 2:
        return u % 8 == 1
    return pow(u % p, (p - 1) // 2, p) == 1


def best_tower(seed, p, max_vmu=6, Lrange=64):
    lev = seed_levels(seed)
    assert len(lev) == 4 and len(set(seed)) == 8
    best = None
    for vmu in range(0, max_vmu + 1):
        mu = p ** vmu
        for L in range(1, Lrange):
            vals = [L + s * mu * r for r in seed for s in (1, -1)]
            if 0 in vals or not all(is_square_Zp(a, p) for a in vals):
                continue
            leaf = max(vp(a, p) // 2 for a in vals)
            ht = max([leaf] + [2 ** (j - 1) * vmu + max(vp(t, p) for t in lev[j - 1])
                               for j in range(1, 5)])
            top = 8 * vmu + vp(lev[3][0], p)
            if best is None or ht < best[0]:
                best = (ht, top, vmu, L)
        if best is not None and best[0] <= 2 ** 3 * vmu:
            break
    return best


if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else 'catalog_1e11.csv'
    primes = [2, 3, 5, 7, 11, 13]
    seeds = seeds_from_catalog(path)
    print('p    min ht_p over seeds (= upper bound on V_p(k,5) for every k)   [top-node val, v(mu), L, seed L1/form]')
    for p in primes:
        cands = []
        for L1, form, s in seeds:
            b = best_tower(s, p)
            if b:
                cands.append((b, L1, form))
        (ht, top, vmu, L), L1, form = min(cands)
        print(f'{p:<4} {ht:<6} [top {top}, v(mu)={vmu}, L={L}, seed {L1} {form}]')
