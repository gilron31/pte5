"""Verify a reconstructed image curve F(J1,J2): at fresh primes, F vanishes at the
J-image of many F_p-points of the family's (s:t:u)-curve (roots via PARI polrootsmod).
Usage: python3 verify_image.py fam1|fam2 FILE
"""
import subprocess
import sys

import sympy as sp

from phi_degree import invariants
from twin_check import CURVES

s, t, u = sp.symbols('s t u')
curve, path = sys.argv[1], sys.argv[2]
L1, L2, L3, _ = invariants(1, t, u)
G = sp.expand(sp.sympify(CURVES[curve]).subs(s, 1))
f = lambda e: str(sp.expand(e)).replace('**', '^')
F = open(path).read().strip().rstrip(',')
gp = f'''
F = {F};
G = {f(G)}; L1 = {f(L1)}; L2 = {f(L2)}; L3 = {f(L3)};
{{
foreach([2147483587, 2147483579, 1000000007], p,
  my(n = 0, bad = 0);
  for (k = 1, 300,
    my(tt = Mod(random(p), p), g = subst(G, 't, tt), R = polrootsmod(lift(g), p));
    for (i = 1, #R,
      my(uu = R[i], a = subst(subst(L1, 't, tt), 'u, uu));
      if (a == 0, next);
      my(j1 = subst(subst(L2, 't, tt), 'u, uu) / a^2, j2 = subst(subst(L3, 't, tt), 'u, uu) / a^4);
      n++; if (subst(subst(F, 'J1, j1), 'J2, j2) != 0, bad++)));
  printf("p = %d: %d curve points, %d where F does not vanish\\n", p, n, bad));
}}
'''
print(subprocess.run(['gp', '-q'], input=gp, capture_output=True, text=True).stdout)
