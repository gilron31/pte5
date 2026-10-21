"""Both Bremner E_4 families, reduced mod p: within- and cross-family twin census.

Family 1 lives on the plane quartic Q4 of Bremner (2-7); there A,B,C,D,E,G are
the shared cubic forms and F,H are recovered from F+H = T and F^2-H^2 = q1-q2.
All eight are multiplied by 2T so that every coordinate is a polynomial.
Family 2 lives on the plane cubic of Bremner Section 3 (explicit forms).
"""
import sys
from collections import defaultdict

import numpy as np
import sympy as sp

from family2 import FORMS as F2_FORMS, CUB, x, y, z

A0 = x**2*y - x*y**2 + 3*x*z**2 - 7*y*z**2
B0 = z*(x**2 + 2*x*y - y**2 + 10*z**2)
C0 = x**2*y - x*y**2 + 7*x*z**2 - 3*y*z**2
D0 = -z*(x**2 - 2*x*y - y**2 - 10*z**2)
E0 = -z*(x**2 + y**2 + 2*x*z + 2*y*z - 10*z**2)
G0 = -z*(x**2 + y**2 - 2*x*z - 2*y*z - 10*z**2)
q1 = (x**4*y**2 - 2*x**3*y**3 + x**2*y**4 + 10*x**3*y*z**2 - 20*x**2*y**2*z**2
      + 10*x*y**3*z**2 - 4*x**3*z**3 - 4*x**2*y*z**3 - 4*x*y**2*z**3 - 4*y**3*z**3
      + 45*x**2*z**4 - 10*x*y*z**4 + 45*y**2*z**4 + 40*x*z**5 + 40*y*z**5)
q2 = q1.subs(z, -z)
T = 2*x**2*y - 2*x*y**2 + 2*x**2*z + 4*x*y*z + 2*y**2*z + 10*x*z**2 - 10*y*z**2
Q4 = (2*x**3*y - 4*x**2*y**2 + 2*x*y**3 + x**3*z + 13*x**2*y*z - 13*x*y**2*z
      - y**3*z + 16*x**2*z**2 + 16*y**2*z**2 + 60*x*z**3 - 60*y*z**3 + 40*z**4)
F1_FORMS = [2*T*A0, 2*T*B0, 2*T*C0, 2*T*D0, 2*T*E0, T**2 + q1 - q2, 2*T*G0, T**2 - q1 + q2]


def check_family1():
    A, B, C, D, E, F, G, H = F1_FORMS
    eqs = [A**2 + B**2 - C**2 - D**2, A**2 + B**2 - E**2 - F**2,
           A**2 + B**2 - G**2 - H**2, A**2*B**2 + C**2*D**2 - E**2*F**2 - G**2*H**2]
    for e in eqs:
        _, r = sp.reduced(sp.expand(e), [Q4], x, y, z, order='grevlex')
        print('  vanishes on Q4:', sp.expand(r) == 0)


def plane_points(poly, p):
    """Projective F_p-points of a plane curve: affine chart z = 1 by vectorized
    Horner in y, then the line at infinity z = 0."""
    coeffs = sp.Poly(poly.subs(z, 1), y).all_coeffs()        # in y, coefficients in x
    cfun = [sp.lambdify(x, c) for c in coeffs]
    ys = np.arange(p, dtype=np.int64)
    for Xv in range(p):
        acc = np.zeros(p, dtype=np.int64)
        for cf in cfun:
            acc = (acc * ys + int(cf(Xv)) % p) % p
        for Yv in np.nonzero(acc == 0)[0]:
            yield Xv, int(Yv), 1
    for Xv, Yv in [(1, t) for t in range(p)] + [(0, 1)]:
        if int(poly.subs({x: Xv, y: Yv, z: 0})) % p == 0:
            yield Xv, Yv, 0


def invariants(vals, p):
    inv = lambda a: pow(a, p - 2, p)
    sq = [int(v) * int(v) % p for v in vals]
    if 0 in sq or len(set(sq)) < 8:
        return None
    L1 = (sq[0] + sq[1]) * inv(2) % p
    if L1 == 0:
        return None
    i8 = inv(8)
    L2 = ((sq[0] - sq[1]) ** 2 + (sq[2] - sq[3]) ** 2) * i8 % p
    u = (sq[0] * sq[1] - sq[2] * sq[3]) % p
    w = (sq[4] * sq[5] - sq[6] * sq[7]) % p
    L3 = (u * u + w * w) * i8 % p
    S = (u * u - w * w) * i8 % p
    l4 = inv(pow(L1, 4, p))
    return (L2 * inv(L1 * L1) % p, L3 * l4 % p), pow(S * l4, 2, p)


def census(p):
    groups = defaultdict(lambda: defaultdict(set))   # (j1,j2) -> family -> {j4}
    counts = {}
    for name, curve, forms in (('F1', Q4, F1_FORMS), ('F2', CUB, F2_FORMS)):
        num = sp.lambdify((x, y, z), forms)
        n = 0
        for P in plane_points(curve, p):
            r = invariants(num(*P), p)
            if r:
                n += 1
                groups[r[0]][name].add(r[1])
        counts[name] = n
    within = {f: sum(1 for g in groups.values() if len(g.get(f, ())) >= 2) for f in ('F1', 'F2')}
    cross = sum(1 for g in groups.values()
                if 'F1' in g and 'F2' in g and (g['F1'] | g['F2']) != (g['F1'] & g['F2']))
    shared = sum(1 for g in groups.values() if 'F1' in g and 'F2' in g)
    return counts, within, shared, cross


if __name__ == '__main__':
    print('symbolic check of family 1 on the quartic:')
    check_family1()
    for p in map(int, sys.argv[1:]):
        counts, within, shared, cross = census(p)
        print(f'p={p}: nondegenerate points {counts}; within-family twins {within}; '
              f'(L1:L2:L3) classes hit by both families {shared}, of which with different L4 {cross}')
