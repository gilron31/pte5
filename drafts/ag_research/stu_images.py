"""Where do family members and the exemplar sit in Bremner's (s:t:u)-plane?"""
import itertools
from fractions import Fraction
from math import gcd

from known_e4 import all_known, twin, primitive
from invert import signed_labelings


def stu_points(r):
    out = set()
    for v in signed_labelings(r):
        A, B, C, D, E, F, G, H = v
        s, t, u = Fraction(C - A, 2), Fraction(A + C, 2), Fraction(E + G, 2)
        out.add(primitive([s, t, u]) if (s, t, u) != (0, 0, 0) else None)
    out.discard(None)
    # normalize sign
    return {p if next(c for c in p if c) > 0 else tuple(-c for c in p) for p in out}


def collinear(P, Q, R):
    M = [P, Q, R]
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
            - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])) == 0


if __name__ == '__main__':
    rows = all_known(kmax=10)
    for fam in ('F1', 'F2'):
        members = [n for n in rows if n.startswith(fam)]
        pts = {n: stu_points(rows[n]) for n in members}
        print(fam, 'images per member:', {n: len(p) for n, p in list(pts.items())[:4]}, '...')
        # lines through at least 4 member images (one image per member, any choice)
        allp = [(n, p) for n, ps in pts.items() for p in ps]
        lines = {}
        for (n1, p1), (n2, p2) in itertools.combinations(allp, 2):
            if n1 == n2:
                continue
            l = (p1[1] * p2[2] - p1[2] * p2[1], p1[2] * p2[0] - p1[0] * p2[2], p1[0] * p2[1] - p1[1] * p2[0])
            g = 0
            for c in l:
                g = gcd(g, c)
            if g == 0:
                continue
            l = tuple(c // g for c in l)
            if next(c for c in l if c) < 0:
                l = tuple(-c for c in l)
            lines.setdefault(l, set()).update({n1, n2})
        best = sorted(lines.items(), key=lambda kv: -len(kv[1]))[:3]
        for l, ns in best:
            print(f'   line {l}: contains images of {len(ns)} members')
    print('exemplar (twin) images:', sorted(stu_points(twin(rows['cat33065996717'])), key=lambda p: max(map(abs, p)))[:6])
