"""Degree of phi: P^2_stu -> P(1,2,4), (s:t:u) -> (L1:L2:L3) on Bremner's slice.

Writes a Singular script computing the fiber of phi through a random point
(s:t:u) = (1:t0:u0), in the chart s = 1, with the L1 = 0 locus saturated away,
and prints the number of fiber points and how they split into Q-components.
Each E_4 on the slice accounts for 8 points of P^2_stu (its signed relabelings);
any further fiber points are different E_4 with the same (L1:L2:L3), i.e. E_5
candidates with both halves on the slice.
Usage: python3 phi_degree.py t0 u0
"""
import subprocess
import sys

import sympy as sp

s, t, u = sp.symbols('s t u')


def invariants(S, T, U):
    A, C, E, G = -S + T, S + T, -S + U, S + U
    b, d = 5*S**2 + 2*S*T + U**2, 5*S**2 - 2*S*T + U**2
    f, h = 5*S**2 + 2*S*U + T**2, 5*S**2 - 2*S*U + T**2
    L1 = A**2 + b                               # 2 L1
    L2 = (A**2 - b)**2 + (C**2 - d)**2          # 8 L2
    W1, W2 = A**2*b - C**2*d, E**2*f - G**2*h
    return [sp.expand(e) for e in (L1, L2, W1**2 + W2**2, W1**2 - W2**2)]  # 8 L3, 8 sqrt(L4)


def script(t0, u0):
    L1, L2, L3, _ = invariants(1, t, u)
    a1, a2, a3, _ = invariants(1, t0, u0)
    e1 = sp.expand(L2 * a1**2 - a2 * L1**2)
    e2 = sp.expand(L3 * a1**4 - a3 * L1**4)
    p = lambda e: str(e).replace('**', '^')
    return f'''LIB "primdec.lib"; LIB "elim.lib";
ring r = 0, (t,u), dp;
poly L1 = {p(L1)};
ideal I = {p(e1)}, {p(e2)};
def S0 = sat(I, ideal(L1));
ideal S; if (typeof(S0) == "list") {{ S = S0[1]; }} else {{ S = S0; }}
S = std(S);
"fiber dim:"; dim(S);
"fiber points with multiplicity:"; vdim(S);
list L = minAssGTZ(S);
"Q-components:"; size(L);
int i;
for (i = 1; i <= size(L); i++) {{ ideal Pi = std(L[i]); print("  deg " + string(vdim(Pi)) + "  mult " + string(vdim(std(L[i]^1)))); }}
"reduced total:"; vdim(std(radical(S)));
quit;
'''


if __name__ == '__main__':
    t0, u0 = int(sys.argv[1]), int(sys.argv[2])
    open('phi_fiber.sing', 'w').write(script(t0, u0))
    out = subprocess.run(['Singular', '-q', 'phi_fiber.sing'], capture_output=True, text=True, timeout=3000)
    print(out.stdout, out.stderr[-2000:])
