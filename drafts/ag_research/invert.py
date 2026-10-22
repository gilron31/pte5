"""Locate a known E_4 on Bremner's slice: recover (x:y:z) and test the curves.

Given tree-ordered roots, try every signed choice (A,B),(C,D),(E,F),(G,H) with
A - C = E - G.  Bremner's chain then gives (lambda:mu:nu:phi) = (s:B:u:D) with
s = (C-A)/2, u = (E+G)/2, and (mu:phi:nu:lambda) is the image of (x:y:z) under
  (10z^2+x^2+2xy-y^2 : 10z^2-x^2+2xy+y^2 : 10z^2-x^2-y^2 : 2z(x+y)).
We solve for (x:y:z) and check that the cubic forms reproduce the roots.
"""
import itertools
import sys
from fractions import Fraction

import sympy as sp

from known_e4 import all_known, twin, primitive
from slice import x, y, z, A_, B_, C_, D_, E_, G_, Q1, Q2, FAM1, FAM2, OTHER

MU = 10*z**2 + x**2 + 2*x*y - y**2
PHI = 10*z**2 - x**2 + 2*x*y + y**2
NU = 10*z**2 - x**2 - y**2
LAM = 2*z*(x + y)


def signed_labelings(r):
    """All (A,B,C,D,E,F,G,H) from root vector r (tree order) with A - C = E - G."""
    pairs = [(r[0], r[1]), (r[2], r[3]), (r[4], r[5]), (r[6], r[7])]
    for half_swap in (False, True):
        P = pairs[2:] + pairs[:2] if half_swap else pairs
        for sw1, sw2 in itertools.product((False, True), repeat=2):
            h1 = [P[1], P[0]] if sw1 else P[:2]
            h2 = [P[3], P[2]] if sw2 else P[2:]
            for flips in itertools.product((0, 1), repeat=4):
                quad = h1 + h2
                lab = []
                for (a, b), f in zip(quad, flips):
                    lab += [b, a] if f else [a, b]
                for signs in itertools.product((1, -1), repeat=8):
                    v = [s * w for s, w in zip(signs, lab)]
                    if v[0] - v[2] == v[4] - v[6] and signs[0] == 1:
                        yield tuple(v)


def locate(v):
    A, B, C, D, E, F, G, H = v
    s, u = Fraction(C - A, 2), Fraction(E + G, 2)
    target = [B, D, u, s]                      # (mu : phi : nu : lambda)
    forms = [MU, PHI, NU, LAM]
    eqs = [sp.expand(forms[i] * target[3] - forms[3] * target[i]) for i in range(3)]
    sols = []
    for zz in (1, 0):
        e = [q.subs(z, zz) for q in eqs]
        try:
            res = sp.solve(e, [x, y], dict=True)
        except Exception:
            res = []
        for so in res:
            X, Y = so.get(x), so.get(y)
            if X is None or Y is None or not (X.is_rational and Y.is_rational):
                continue
            sols.append((sp.Rational(X), sp.Rational(Y), zz))
    out = []
    for X, Y, Z in sols:
        den = sp.ilcm(X.q, Y.q)
        P = (int(X * den), int(Y * den), int(Z * den))
        vals = [int(f.subs({x: P[0], y: P[1], z: P[2]})) for f in (A_, B_, C_, D_, E_, G_)]
        if all(vals[i] * A == vals[0] * w for i, w in enumerate((A, B, C, D, E, G))) and vals[0]:
            out.append(P)
    return out


def curves_at(P):
    d = {'FAM1': FAM1, 'FAM2': FAM2, **OTHER}
    return {n: int(f.subs({x: P[0], y: P[1], z: P[2]})) for n, f in d.items()}


if __name__ == '__main__':
    rows = all_known(kmax=3)
    names = sys.argv[1:] or ['cat33065996717', 'cat43937864165', 'cat67405', 'cat5915065']
    for name in names:
        for form, r in (('root', rows[name]), ('twin', twin(rows[name]))):
            pts = set()
            for v in signed_labelings(r):
                pts.update(locate(v))
            for P in sorted(pts):
                g = sp.igcd(sp.igcd(*P[:2]), P[2])
                P = tuple(c // g for c in P)
                vals = curves_at(P)
                hit = [n for n, val in vals.items() if val == 0]
                print(f'{name} ({form} form): (x:y:z) = {P}  on {hit or "NO KNOWN CURVE"}')
