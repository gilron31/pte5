"""Preimages, on a family's (s:t:u) curve, of the rational singular points of
its (J1, J2) image; for each, are they rational slice points (b, d, f, h
squares) giving different E_4 (different L4) with a square L1 ratio?
A within-family twin (an E_5 from two members of one family) must show up here.
Usage: python3 twin_check.py fam1|fam2 'J1a J2a; J1b J2b; ...'
"""
import subprocess
import sys
from fractions import Fraction

import sympy as sp

from phi_degree import invariants

s, t, u = sp.symbols('s t u')
CURVES = {
    'fam1': "-500*s**4*t**2+61*s**2*t**4+16*t**6+1500*s**4*t*u-436*s**2*t**3*u-32*t**5*u-500*s**4*u**2+806*s**2*t**2*u**2-92*t**4*u**2-436*s**2*t*u**3+236*t**3*u**3+61*s**2*u**4-92*t**2*u**4-32*t*u**5+16*u**6",
    'fam2': "10000*s**7+500*s**5*t**2+500*s**4*t**3-400*s**3*t**4-25*s*t**6-3000*s**5*t*u-100*s**4*t**2*u-200*s**3*t**3*u-40*s**2*t**4*u+70*s*t**5*u+500*s**5*u**2-100*s**4*t*u**2+800*s**3*t**2*u**2-120*s**2*t**3*u**2+29*s*t**4*u**2-12*t**5*u**2+500*s**4*u**3-200*s**3*t*u**3-120*s**2*t**2*u**3-68*s*t**3*u**3+28*t**4*u**3-400*s**3*u**4-40*s**2*t*u**4+29*s*t**2*u**4+28*t**3*u**4+70*s*t*u**5-12*t**2*u**5-25*s*u**6",
}


def is_sq(q):
    q = Fraction(q)
    if q < 0:
        return False
    n, d = q.numerator, q.denominator
    from math import isqrt
    return isqrt(n) ** 2 == n and isqrt(d) ** 2 == d


def preimages(curve, J1, J2):
    G = sp.sympify(CURVES[curve]).subs(s, 1)
    L1, L2, L3, _ = invariants(1, t, u)
    def p(e):
        poly = sp.Poly(sp.expand(e), t, u)
        _, poly = poly.clear_denoms()
        return str(poly.as_expr()).replace('**', '^')
    src = f'''LIB "primdec.lib";
ring r = 0, (t,u), dp;
poly L1 = {p(L1)};
ideal I = {p(G)}, {p(sp.Rational(J1) * L1**2 - L2)}, {p(sp.Rational(J2) * L1**4 - L3)};
def S0 = sat(I, ideal(L1));
ideal S; if (typeof(S0) == "list") {{ S = S0[1]; }} else {{ S = S0; }}
S = std(S);
"LEN"; vdim(S);
list L = minAssGTZ(S);
int i; for (i = 1; i <= size(L); i++) {{ ideal P = std(L[i]); if (vdim(P) == 1) {{ "PT"; P[1]; P[2]; }} else {{ "CLUSTER"; vdim(P); }} }}
quit;
'''
    open('twin_tmp.sing', 'w').write(src)
    out = subprocess.run(['Singular', '-q', 'twin_tmp.sing'], capture_output=True, text=True).stdout
    lines = [l.strip() for l in out.splitlines() if 'redefin' not in l]
    pts, clusters, length = [], [], None
    i = 0
    while i < len(lines):
        if lines[i] == 'LEN':
            length = int(lines[i + 1]); i += 2
        elif lines[i] == 'PT':
            eqs = [sp.sympify(__import__('re').sub(r'(\d)([tu])', r'\1*\2', lines[i + k].replace('^', '**'))) for k in (1, 2)]
            so = sp.solve(eqs, [t, u], dict=True)[0]
            pts.append((Fraction(str(so[t])), Fraction(str(so[u])))); i += 3
        elif lines[i] == 'CLUSTER':
            clusters.append(int(lines[i + 1])); i += 2
        else:
            i += 1
    return length, pts, clusters


if __name__ == '__main__':
    curve = sys.argv[1]
    for item in sys.argv[2].split(';'):
        J1, J2 = item.split()
        length, pts, clusters = preimages(curve, J1, J2)
        print(f'(J1, J2) = ({J1}, {J2}): fiber length {length}; rational points {len(pts)}; clusters {clusters}')
        for T, U in pts:
            S, T_, U_ = Fraction(1), T, U
            b, d = 5 + 2*T_ + U_**2, 5 - 2*T_ + U_**2
            f, h = 5 + 2*U_ + T_**2, 5 - 2*U_ + T_**2
            L1, L2, L3, S4 = invariants(S, T_, U_)
            on_slice = all(is_sq(v) for v in (b, d, f, h))
            j4 = Fraction(S4) ** 2 / Fraction(L1) ** 8
            print(f'    (1:{T}:{U})  slice point: {on_slice}  j4 = {j4}  L1 = {L1}')
