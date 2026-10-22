"""Q1 ansatz: family-2-type curves through a given slice point P0.

Family 2 satisfies q1 - F2^2 = (y+z) * conic * FAM2 and q2 - H2^2 = 4z(x-z)(y+z) * FAM2,
with F2 = Fy + (y+z) Q, H2 = Hy + (y+z) Q' for small quadrics Q, Q' (Fy, Hy the
square roots of q1, q2 on the split line y = -z). Search: integral quadrics Q with
Phi = e2*Fy + (y+z) Q and Phi(P0) = e1*sqrt(q1(P0)), coefficients at most B, via
PARI qfminim on a Kannan embedding; factor (q1 - Phi^2)/(y+z) in Singular; keep the
cubic factors through P0 on which q2 is a square at every F_p-point (several p).
Usage: python3 ansatz_search.py x0 y0 z0 [B]
"""
import itertools
import re
import subprocess
import sys
from math import isqrt

import sympy as sp

from slice import Q1, Q2, x, y, z

FY = z*(x**2 - 6*x*z - 3*z**2)
HY = z*(x**2 - 2*x*z + 9*z**2)
# split lines l = 0 with q1|_l = R1^2, q2|_l = R2^2 (R written in the remaining variables)
LINES = {
    'y+z': (y + z, z*(x**2 - 6*x*z - 3*z**2), z*(x**2 - 2*x*z + 9*z**2)),
    'y-z': (y - z, z*(x**2 + 2*x*z + 9*z**2), z*(x**2 + 6*x*z - 3*z**2)),
    'x+z': (x + z, z*(y**2 - 6*y*z - 3*z**2), z*(y**2 - 2*y*z + 9*z**2)),
    'x-z': (x - z, z*(y**2 + 2*y*z + 9*z**2), z*(y**2 + 6*y*z - 3*z**2)),
    'x+y': (x + y, 2*y*(y**2 - 5*z**2), 2*y*(y**2 - 5*z**2)),
    'z': (z, x*y*(x - y), x*y*(x - y)),
}
MONS = [x**2, y**2, z**2, x*y, x*z, y*z]


def lattice_candidates(P0, target, B):
    """Integral c with sum c_i m_i(P0) = target and max|c_i| <= B (via PARI)."""
    m = [int(mon.subs({x: P0[0], y: P0[1], z: P0[2]})) for mon in MONS]
    gp = f'''
m = {m}; T = {target}; B = {B};
g = content(m);
if (T % g, print("NONE"); quit);
c0 = matsolvemod(Mat(m), [abs(T)]~ * 0 + T, [T]~ * 0 + 0, 1);
'''
    # simpler and robust: particular solution by extended gcd chain, kernel by matkerint
    gp = f'''
m = {m}; T = {target}; B = {B};
{{
my(g = content(m));
if (T % g != 0, print("NONE"); quit);
my(K = matkerint(Mat(m)), c0, v = vector(#m), r = T, M, A, L);
\\\\ particular solution: c0 = lift of T along a Bezout combination
my(G = matsnf(Mat(m), 1));
c0 = matsolve(Mat(m)~ * Mat(m), Mat(m)~ * T);
}}
'''
    # use a cleaner approach in PARI: solve via the HNF of the 1 x 6 matrix
    gp = f'''
m = {m}; T = {target}; B = {B};
{{
my(g = content(m), K, U, c0, M, A, S, n = 0);
if (T % g != 0, print("NONE"); quit);
my(H = mathnf(Mat(m), 1));       \\\\ [H, U] with Mat(m) * U = [0,...,0, g]
U = H[2];
c0 = U[, #m] * (T / g);
K = matrix(#m, #m - 1, i, j, U[i, j]);
K = K * qflll(K);
M = B;
A = matconcat([K, c0; matrix(1, #m - 1), Mat(M)]);
S = qfminim(A~ * A, 6 * B^2 + M^2, 200000, 2)[3];
for (i = 1, #S[1, ],
  my(w = A * S[, i]);
  if (abs(w[#w]) != M, next);
  my(c = vector(#m, j, w[j] * sign(w[#w])));
  if (vecmax(apply(abs, c)) <= B, n++; print(c)));
print("COUNT ", n);
}}
'''
    out = subprocess.run(['gp', '-q'], input=gp, capture_output=True, text=True, timeout=600).stdout
    cands = []
    for line in out.splitlines():
        line = line.strip()
        if line.startswith('['):
            cands.append([int(v) for v in line.strip('[]').split(',')])
    return cands


def square_on_curve(C, q, primes=(101, 103, 107, 109, 113)):
    """Is q a square at every F_p point of the plane curve C (several p)? (necessary test)"""
    Cp = sp.Poly(C, x, y, z)
    for p in primes:
        cnt = 0
        for X in range(p):
            for Y in range(p):
                for Z in ((1,) if True else ()):
                    if Cp.eval((X, Y, Z)) % p == 0:
                        v = int(q.subs({x: X, y: Y, z: Z})) % p
                        if v and pow(v, (p - 1) // 2, p) != 1:
                            return False
                        cnt += 1
        if cnt == 0:
            continue
    return True


def run(P0, B=12, line="y+z"):
    ell, R1, _ = LINES[line]
    q1v = int(Q1.subs({x: P0[0], y: P0[1], z: P0[2]}))
    q2v = int(Q2.subs({x: P0[0], y: P0[1], z: P0[2]}))
    F0, H0 = isqrt(q1v), isqrt(q2v)
    assert F0 * F0 == q1v and H0 * H0 == q2v, 'not a slice point'
    w = int(ell.subs({x: P0[0], y: P0[1], z: P0[2]}))
    fy = int(R1.subs({x: P0[0], y: P0[1], z: P0[2]}))
    found = []
    quintics = []
    for e1, e2 in itertools.product((1, -1), repeat=2):
        num = e1 * F0 - e2 * fy
        if w == 0 or num % w:
            continue
        for c in lattice_candidates(P0, num // w, B):
            Q = sum(ci * mi for ci, mi in zip(c, MONS))
            Phi = e2 * R1 + ell * Q
            R = sp.cancel((Q1 - Phi**2) / ell)
            quintics.append((e1, e2, tuple(c), sp.expand(R)))
    if not quintics:
        return []
    # factor all quintics in one Singular run
    p = lambda e: str(e).replace('**', '^')
    lines = ['ring r = 0, (x,y,z), dp;', 'short = 0;', 'list L; poly R;']
    for i, (_, _, _, R) in enumerate(quintics):
        lines.append(f'R = {p(R)}; L = factorize(R); print("F{i}"); int j; for (j = 1; j <= size(L[1]); j++) {{ if (deg(L[1][j]) == 3) {{ print(L[1][j]); }} }} kill j;')
    lines.append('quit;')
    open('ansatz_tmp.sing', 'w').write('\n'.join(lines))
    out = subprocess.run(['Singular', '-q', 'ansatz_tmp.sing'], capture_output=True, text=True, timeout=1800).stdout
    cur = None
    cubics = {}
    for line in out.splitlines():
        line = line.strip()
        if re.fullmatch(r'F\d+', line):
            cur = int(line[1:])
        elif line and cur is not None and any(v in line for v in 'xyz'):
            expr = sp.sympify(re.sub(r'(\d)([xyz])', r'\1*\2', line).replace('^', '**'))
            if expr.subs({x: P0[0], y: P0[1], z: P0[2]}) == 0:
                cubics.setdefault(sp.factor(expr), []).append(quintics[cur][:3])
    for C, src in cubics.items():
        ok = square_on_curve(C, Q2, primes=(101, 103))
        found.append((C, ok, src[:2]))
    return len(quintics), found


if __name__ == '__main__':
    P0 = tuple(int(v) for v in sys.argv[1:4])
    B = int(sys.argv[4]) if len(sys.argv) > 4 else 12
    line = sys.argv[5] if len(sys.argv) > 5 else 'y+z'
    res = run(P0, B, line)
    if not res:
        print(P0, line, 'no admissible sign/divisibility choice')
    else:
        n, found = res
        print(f'{P0} [{line}]: {n} candidate quadrics; cubic factors through P0: {len(found)}')
        for C, ok, src in found:
            print(f'   {C}   q2 square on it mod 101, 103: {ok}   from {src}')
