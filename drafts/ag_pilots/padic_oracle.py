"""An independent tree-walk oracle for towers over Z/p^k, and explicit p-adic E_5s.

Construction (Corollary 3.2 read over Z_p): take an E_4 seed with roots rho
(tree order) and a scale mu with v_p(mu) >= 1 (p odd).  Then L = 1 makes every
1 +- mu*rho_k = 1 mod p a square in Z_p, so
    (L_1..L_5) = (1, mu^2 L1', mu^4 L2', mu^8 L3', mu^16 L4')
is a nondegenerate E_5 over Z_p.  Its node valuations are bounded, so the
valuation floor V_p(k,5) of paper.md section 4 is bounded for every k.

The oracle walks the tree of nested radicals one tuple at a time (the
definition, no sets): splits(j, c) <=> exists s, s^2 = L_j + c (mod p^k),
v(s) <= cap, splits(j-1, s) and splits(j-1, -s); splits(0, c) = True.
"""
import csv
import sys
from fractions import Fraction
from functools import lru_cache
from math import isqrt


def params_from_roots(roots):
    """(L_1..L_d) of the tower with these tree-ordered positive roots (rationals)."""
    out = []
    level = [Fraction(r) for r in roots]
    while True:
        if len(level) == 1:
            out.append(level[0] ** 2)
            return out
        out.append((level[0] ** 2 + level[1] ** 2) / 2)
        level = [(level[2 * k] ** 2 - level[2 * k + 1] ** 2) / 2 for k in range(len(level) // 2)]


def check_tree(roots):
    """All tree conditions hold exactly (each level's pairs have equal sum of squares)."""
    level = [Fraction(r) for r in roots]
    while len(level) > 1:
        sums = {level[2 * k] ** 2 + level[2 * k + 1] ** 2 for k in range(len(level) // 2)}
        if len(sums) != 1:
            return False
        level = [(level[2 * k] ** 2 - level[2 * k + 1] ** 2) / 2 for k in range(len(level) // 2)]
    return True


class Oracle:
    def __init__(self, p, k):
        self.p, self.k, self.m = p, k, p ** k
        self.roots = {}
        for x in range(self.m):
            self.roots.setdefault(x * x % self.m, []).append(x)
        self.val = [k] + [self._v(x) for x in range(1, self.m)]

    def _v(self, x):
        v = 0
        while x % self.p == 0:
            x //= self.p
            v += 1
        return v

    def red(self, q):
        q = Fraction(q)
        return q.numerator * pow(q.denominator, -1, self.m) % self.m

    def floor(self, ells):
        L = [self.red(x) for x in ells]
        d = len(L)
        for cap in range(self.k + 1):
            @lru_cache(maxsize=None)
            def splits(j, c):
                if j == 0:
                    return True
                for s in self.roots.get((L[j - 1] + c) % self.m, ()):
                    if self.val[s] <= cap and splits(j - 1, s) and splits(j - 1, (-s) % self.m):
                        return True
                return False
            if splits(d, 0):
                return cap
        return self.k


EXEMPLAR = [33065996717, 489628056848329146064,
            175480010455650701584492675662518592000000,
            40042900368028062136207226327668992000000 ** 2]


def catalog_point_seeds(path):
    for row in csv.DictReader(open(path)):
        L1 = int(row['l1'])
        pts = []
        for i in range(4):
            nu = int(row[f'nu{i}'])
            pts.append((isqrt((L1 + nu) // 2), isqrt((L1 - nu) // 2)))
        yield L1, [v for m, n in pts for v in (m, n)]


if __name__ == '__main__':
    cat = sys.argv[1] if len(sys.argv) > 1 else 'catalog_1e11.csv'
    print('calibration: the paper.md section 7 exemplar E_4 (paper: ht_3 = 6, ht_7 = 2)')
    for p, k in ((3, 7), (3, 8), (7, 3), (7, 4)):
        print(f'  p={p} k={k}: tuple floor {Oracle(p, k).floor(EXEMPLAR)}')
    seeds = dict(catalog_point_seeds(cat))
    seed = seeds[67405]
    assert check_tree(seed)
    ells = [Fraction(1)] + params_from_roots([3 * r for r in seed])
    print('E_5 over Z_3 from seed 67405 (point form), mu = 3, L = 1: ht_3 = 14, so below k = 15')
    print('the tuple reads as degenerate (floor = k) -- the precision effect of paper.md 6.6:')
    for k in (8, 9, 10):
        print(f'  p=3 k={k}: tuple floor {Oracle(3, k).floor(ells)}')
