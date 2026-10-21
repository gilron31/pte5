"""Bremner's second E_4 family: symbolic check, then a mod-p twin census.

A "twin" is a pair of family members with the same (L1:L2:L3) in weighted
projective space but different L4 -- i.e. an E_5 made of two family members.
If the map pi: C -> P(1,2,4) had a twin component, mod-p collisions would grow
like p; isolated collisions (nodes of the image curve) stay bounded.
"""
import sys
from collections import defaultdict

import sympy as sp

x, y, z = sp.symbols('x y z')
CUB = 3*x**2*y - 4*x*y**2 + y**3 + (x**2 - 14*x*y + 9*y**2)*z + 10*(x + y)*z**2
FORMS = [
    x**2*y - x*y**2 + (3*x - 7*y)*z**2,                                   # A
    z*(x**2 + 2*x*y - y**2 + 10*z**2),                                    # B
    x**2*y - x*y**2 + (7*x - 3*y)*z**2,                                   # C
    -z*(x**2 - 2*x*y - y**2 - 10*z**2),                                   # D
    -z*(x**2 + y**2 + 2*(x + y)*z - 10*z**2),                             # E
    -2*x**2*y + 3*x*y**2 - y**3 - (x**2 - 8*x*y + 7*y**2)*z - (x + 3*y)*z**2,  # F
    -z*(x**2 + y**2 - 2*(x + y)*z - 10*z**2),                             # G
    -x**2*y + x*y**2 + (6*x*y - 2*y**2)*z + (3*x - 11*y)*z**2,            # H
]


def symbolic_check():
    A, B, C, D, E, F, G, H = FORMS
    eqs = [A**2 + B**2 - C**2 - D**2, A**2 + B**2 - E**2 - F**2,
           A**2 + B**2 - G**2 - H**2, A**2*B**2 + C**2*D**2 - E**2*F**2 - G**2*H**2]
    for e in eqs:
        _, r = sp.reduced(sp.expand(e), [CUB], x, y, z, order='grevlex')
        print('  vanishes on the cubic:', sp.expand(r) == 0)


NUM = sp.lambdify((x, y, z), FORMS)


def points(p):
    """Projective points of the cubic over F_p (it is quadratic in z)."""
    inv = lambda a: pow(a, p - 2, p)
    for X, Y in [(1, t) for t in range(p)] + [(0, 1)]:
        a = 10 * (X + Y) % p
        b = (X * X - 14 * X * Y + 9 * Y * Y) % p
        c = (3 * X * X * Y - 4 * X * Y * Y + Y ** 3) % p
        if a == 0:
            if b:
                yield X, Y, (-c * inv(b)) % p
            continue
        disc = (b * b - 4 * a * c) % p
        if disc == 0:
            yield X, Y, (-b * inv(2 * a)) % p
        elif pow(disc, (p - 1) // 2, p) == 1:
            rt = sp.sqrt_mod(disc, p)
            for s in (rt, p - rt):
                yield X, Y, ((-b + s) * inv(2 * a)) % p
    if 0 == 0:  # the point (0:0:1) satisfies the cubic
        yield 0, 0, 1


def census(p):
    inv = lambda a: pow(a, p - 2, p)
    groups = defaultdict(set)
    npts = 0
    for X, Y, Z in points(p):
        npts += 1
        sq = [int(v) * int(v) % p for v in NUM(X, Y, Z)]
        if 0 in sq or len(set(sq)) < 8:
            continue  # degenerate mod p
        L1 = (sq[0] + sq[1]) * inv(2) % p
        if L1 == 0:
            continue
        i8 = inv(8)
        L2 = ((sq[0] - sq[1]) ** 2 + (sq[2] - sq[3]) ** 2) * i8 % p
        u = (sq[0] * sq[1] - sq[2] * sq[3]) % p
        w = (sq[4] * sq[5] - sq[6] * sq[7]) % p
        L3 = (u * u + w * w) * i8 % p
        S = (u * u - w * w) * i8 % p
        l4 = inv(pow(L1, 4, p))
        groups[(L2 * inv(L1 * L1) % p, L3 * l4 % p)].add(pow(S * l4, 2, p))
    twins = sum(1 for s in groups.values() if len(s) >= 2)
    return npts, len(groups), twins


if __name__ == '__main__':
    print('symbolic check of family 2 on the cubic:')
    symbolic_check()
    for p in map(int, sys.argv[1:]):
        n, g, t = census(p)
        print(f'p={p}: points {n}, distinct (L1:L2:L3) {g}, classes with >=2 distinct L4: {t}')
