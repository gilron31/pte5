"""Bremner's slice S_B = V_4 cap {A - C = E - G}, in his (x:y:z) coordinates.

A point (x:y:z) of P^2 gives A,B,C,D,E,G as cubic forms; it is an E_4 exactly
when q1 = F^2 and q2 = H^2 are both squares (Bremner 2008, eqs. 2-5, 2-6).
The two known families are the curves FAM1 (quartic, rank 1) and FAM2 (cubic,
rank 1); Bremner's elimination also produces the factors listed in OTHER.
"""
from math import gcd, isqrt

import sympy as sp

x, y, z = sp.symbols('x y z')

A_ = x**2*y - x*y**2 + 3*x*z**2 - 7*y*z**2
B_ = z*(x**2 + 2*x*y - y**2 + 10*z**2)
C_ = x**2*y - x*y**2 + 7*x*z**2 - 3*y*z**2
D_ = -z*(x**2 - 2*x*y - y**2 - 10*z**2)
E_ = -z*(x**2 + y**2 + 2*x*z + 2*y*z - 10*z**2)
G_ = -z*(x**2 + y**2 - 2*x*z - 2*y*z - 10*z**2)
Q1 = (x**4*y**2 - 2*x**3*y**3 + x**2*y**4 + 10*x**3*y*z**2 - 20*x**2*y**2*z**2
      + 10*x*y**3*z**2 - 4*x**3*z**3 - 4*x**2*y*z**3 - 4*x*y**2*z**3 - 4*y**3*z**3
      + 45*x**2*z**4 - 10*x*y*z**4 + 45*y**2*z**4 + 40*x*z**5 + 40*y*z**5)
Q2 = Q1.subs(z, -z)

FAM1 = (2*x**3*y - 4*x**2*y**2 + 2*x*y**3 + x**3*z + 13*x**2*y*z - 13*x*y**2*z
        - y**3*z + 16*x**2*z**2 + 16*y**2*z**2 + 60*x*z**3 - 60*y*z**3 + 40*z**4)
FAM2 = 3*x**2*y - 4*x*y**2 + y**3 + (x**2 - 14*x*y + 9*y**2)*z + 10*(x + y)*z**2
OTHER = {
    'rank0_cubic': x**2*y - x*y**2 + x**2*z + y**2*z + 5*x*z**2 - 5*y*z**2 - 10*z**3,
    'x+y': x + y, 'x-z': x - z, 'z': z, 'y+z': y + z,
}

_forms = sp.lambdify((x, y, z), [A_, B_, C_, D_, E_, G_, Q1, Q2], 'math')
_curves = {name: sp.lambdify((x, y, z), f, 'math')
           for name, f in [('FAM1', FAM1), ('FAM2', FAM2)] + list(OTHER.items())}


def e4_at(X, Y, Z):
    """The E_4 (A..H) at an integer point, or None if q1 or q2 is not a square."""
    A, B, C, D, E, G, q1, q2 = (int(v) for v in _forms(X, Y, Z))
    if q1 < 0 or q2 < 0:
        return None
    F, H = isqrt(q1), isqrt(q2)
    if F * F != q1 or H * H != q2:
        return None
    return [A, B, C, D, E, F, G, H]


def canonical(roots):
    """Primitive sorted |roots|, or None if degenerate (a zero or a repeated |root|)."""
    r = [abs(v) for v in roots]
    if 0 in r or len(set(r)) < 8:
        return None
    g = 0
    for v in r:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in r))


def on_curves(X, Y, Z):
    return [name for name, f in _curves.items() if f(X, Y, Z) == 0]


def tower(roots):
    """(L1, L2, L3, sqrt(L4)) of the E_4 with these roots in Bremner's pairing
    (A,B),(C,D) | (E,F),(G,H), as exact rationals."""
    from fractions import Fraction
    A, B, C, D, E, F, G, H = (Fraction(v) for v in roots)
    L1 = (A**2 + B**2) / 2
    L2 = ((A**2 - B**2)**2 + (C**2 - D**2)**2) / 8
    u, w = A**2 * B**2 - C**2 * D**2, E**2 * F**2 - G**2 * H**2
    return L1, L2, (u * u + w * w) / 8, (u * u - w * w) / 8
