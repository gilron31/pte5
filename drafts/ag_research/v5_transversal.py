"""B1 (fast): transversal type of V_5 along each curve stratum of N12.

For each K with a 1-dimensional singular stratum, over F_p (p = 32003):
  1. stratum ideal in the chart sigma = 1 (v5_singular.py's criterion), plus a
     random affine hyperplane -> finitely many points; keep an F_p-rational one;
  2. local model X_K = { V_4(d) = 0, 2 y_k^2 = 1 - d_k (k in K) } at that point
     (y_K = 0), cut by a random hyperplane through it (transversal to the curve);
  3. length of the cut's singular scheme at the point (local ordering).
Length 1 means a transversal ordinary double point, i.e. cA_1 along the curve.
Usage: python3 v5_transversal.py
"""
import random
import re
import subprocess

from v5_singular import V4

P = 32003
CURVES = [(0, 2), (0, 4), (0, 1, 2, 3), (0, 1, 4, 5), (0, 2, 4, 6)]


def run(src):
    return subprocess.run(['Singular', '-q'], input=src, capture_output=True, text=True, timeout=600).stdout


def stratum_points(K, rng, tries=12):
    cols = [k for k in range(8) if k not in K]
    for _ in range(tries):
        lin = ' + '.join(f'{rng.randrange(1, P)}*d{k}' for k in range(8)) + f' - {rng.randrange(P)}'
        src = f'''LIB "primdec.lib"; LIB "elim.lib";
ring r = {P}, (d0,d1,d2,d3,d4,d5,d6,d7), dp;
short = 0;
ideal V = {", ".join(V4).replace("sg", "1")};
matrix N = jacob(V);
matrix M[4][{len(cols) + 1}];
int i;
''' + '\n'.join(f'for (i = 1; i <= 4; i++) {{ M[i, {c + 1}] = N[i, {k + 1}]; }}' for c, k in enumerate(cols)) + f'''
for (i = 1; i <= 4; i++) {{ M[i, {len(cols) + 1}] = {" + ".join(f"N[i, {k + 1}]" for k in K)}; }}
ideal I = V{"".join(f", d{k} - 1" for k in K)}, minor(M, 4), {lin};
poly opn = {"*".join(f"(1 - d{k}^2)" for k in cols)};
def S0 = sat(I, ideal(opn));
ideal S; if (typeof(S0) == "list") {{ S = S0[1]; }} else {{ S = S0; }}
list L = minAssGTZ(S);
int j;
for (j = 1; j <= size(L); j++) {{ ideal Q = std(L[j]); if (vdim(Q) == 1) {{ "PT"; print(Q); "END"; }} }}
quit;
'''
        out = run(src)
        m = re.search(r'PT\n(.*?)END', out, re.S)
        if m:
            pt = {}
            for line in m.group(1).replace(',', '\n').splitlines():
                line = line.strip()
                mm = re.fullmatch(r'(?:_\[\d+\]=)?(d\d)([+-]\d+)?', line)
                if mm:
                    pt[mm.group(1)] = (-int(mm.group(2) or 0)) % P
            if len(pt) == 8:
                return [pt[f'd{k}'] for k in range(8)]
    return None


def transversal_length(K, pt, rng):
    ys = [f'y{k}' for k in K]
    shift = [f'poly d{k} = {pt[k]} + e{k};' for k in range(8)]
    eqs = [v.replace('sg', '1') for v in V4] + [f'2*y{k}^2 - (1 - d{k})' for k in K]
    cut = ' + '.join(f'{rng.randrange(1, P)}*e{k}' for k in range(8)) + ''.join(
        f' + {rng.randrange(1, P)}*y{k}' for k in K)
    n = len(eqs) + 1
    src = f'''ring r = {P}, ({",".join(ys)},e0,e1,e2,e3,e4,e5,e6,e7), ds;
{" ".join(shift)}
ideal I = {", ".join(eqs)}, {cut};
"ON"; size(jet(I, 0));
ideal S = std(I + minor(jacob(I), {n}));
"LEN"; vdim(S);
quit;
'''
    out = run(src)
    on = re.search(r'ON\n(\d+)', out)
    ln = re.search(r'LEN\n(-?\d+)', out)
    return (on and on.group(1) == '0'), int(ln.group(1)) if ln else None


if __name__ == '__main__':
    rng = random.Random(7)
    for K in CURVES:
        pt = stratum_points(K, rng)
        if pt is None:
            print(f'K = {K}: no F_{P}-point found on the stratum'); continue
        on, ln = transversal_length(K, pt, rng)
        verdict = 'transversal ODP (cA1, canonical)' if ln == 1 else f'length {ln}: NOT an ODP, look closer'
        print(f'K = {K!s:<14} point {pt}: on X: {on}; {verdict}', flush=True)
