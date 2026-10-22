"""The rational points of a phi-fiber, and whether they are the same E_4."""
import sys
import subprocess
from fractions import Fraction
import sympy as sp
from phi_degree import invariants, script

t0, u0 = int(sys.argv[1]), int(sys.argv[2])
src = script(t0, u0).replace('list L = minAssGTZ(S);', 'list L = minAssGTZ(S);\nint j; for (j = 1; j <= size(L); j++) { if (vdim(std(L[j])) == 1) { print(std(L[j])); } }')
open('phi_fiber2.sing', 'w').write(src)
out = subprocess.run(['Singular', '-q', 'phi_fiber2.sing'], capture_output=True, text=True).stdout
t, u = sp.symbols('t u')
pts = []
for line in out.splitlines():
    pass
# parse the printed linear ideals "_[1]=u-..., _[2]=t-..."
cur = []
import re
for line in out.splitlines():
    line = line.strip().rstrip(',')
    if re.fullmatch(r'[-+0-9tu/ ]+', line) and ('t' in line or 'u' in line):
        cur.append(sp.sympify(re.sub(r'(\d)([tu])', r'\1*\2', line)))
        if len(cur) == 2:
            so = sp.solve(cur, [t, u], dict=True)[0]
            pts.append((Fraction(str(so[t])), Fraction(str(so[u]))))
            cur = []
print(len(pts), 'rational fiber points')
def sqroots(S, T, U):
    A, C, E, G = -S + T, S + T, -S + U, S + U
    b, d = 5*S**2 + 2*S*T + U**2, 5*S**2 - 2*S*T + U**2
    f, h = 5*S**2 + 2*S*U + T**2, 5*S**2 - 2*S*U + T**2
    sq = [A*A, b, C*C, d, E*E, f, G*G, h]
    L1 = (A*A + b) / 2
    return tuple(sorted(x / L1 for x in sq)), invariants(S, T, U)
groups = {}
for T, U in pts:
    key, inv = sqroots(Fraction(1), T, U)
    L1, L2, L3, S4 = inv
    j4 = S4 * S4 / L1**8
    groups.setdefault(key, []).append((T, U, j4))
print(len(groups), 'distinct E_4 (normalized squared-root multisets) among them')
for k, v in groups.items():
    print('  ', len(v), 'points; j4 values', sorted({str(x[2])[:30] for x in v}))
