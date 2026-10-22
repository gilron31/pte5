"""G1.2: the singular locus of V_5, stratum by stratum (Kummer-tower criterion).

V_5 lies over the cone Y = {(sigma : d) : d in V_4} in P^8 via
s_{2k} = (sigma + d_k)/2, s_{2k+1} = (sigma - d_k)/2, r_j^2 = s_j.
Chain rule through r -> s: at a point with sigma != 0 whose zero leaves are
r_{2k+1} = 0 for k in K (so d_k = sigma there; signs are a symmetry), V_5 is
singular iff the 4 x (9 - |K|) matrix
      N'_K(d) = [ columns of N(d) outside K  |  sum of the columns of N(d) in K ]
has rank < 4, where N(d) is the 4 x 8 Jacobian of V_4's equations in d.
(K empty: rank N(d) < 4, i.e. the cone over Sing(V_4).)  sigma = 0 strata
(L_1 = 0) are degenerate and non-real and are skipped.

For each K up to the tree symmetries, Singular computes the stratum
{d in V_4, d_k = sigma (k in K), d_k != +-sigma (k not in K), sigma != 0, rank < 4}.
Usage: python3 v5_singular.py
"""
import itertools
import subprocess

# tree symmetries of the 8 leaves: swap within a pair, swap pairs in a half, swap halves
def tree_group():
    gens = [(1, 0, 2, 3, 4, 5, 6, 7), (2, 3, 0, 1, 4, 5, 6, 7), (4, 5, 6, 7, 0, 1, 2, 3)]
    G = {tuple(range(8))}
    frontier = list(G)
    while frontier:
        g = frontier.pop()
        for h in gens:
            gh = tuple(g[h[i]] for i in range(8))
            if gh not in G:
                G.add(gh)
                frontier.append(gh)
    return G


def orbit_reps():
    G = tree_group()
    seen, reps = set(), []
    for n in range(0, 9):
        for K in itertools.combinations(range(8), n):
            key = min(tuple(sorted(g[k] for k in K)) for g in G)
            if key not in seen:
                seen.add(key)
                reps.append(key)
    return reps, len(G)


V4 = ['d0^2+d1^2-d2^2-d3^2', 'd0^2+d1^2-d4^2-d5^2', 'd0^2+d1^2-d6^2-d7^2',
      'd0^2*d1^2+d2^2*d3^2-d4^2*d5^2-d6^2*d7^2']


def stratum_script(K):
    cols = [k for k in range(8) if k not in K]
    lines = ['LIB "primdec.lib"; LIB "elim.lib";',
             'ring r = 0, (sg,d0,d1,d2,d3,d4,d5,d6,d7), dp;',
             'ideal V = ' + ', '.join(V4) + ';',
             'matrix N = jacob(V);',          # 4 x 9 (column 1 is d/dsigma = 0)
             f'matrix M[4][{len(cols) + (1 if K else 0)}];',
             'int i;']
    for c, k in enumerate(cols):
        lines.append(f'for (i = 1; i <= 4; i++) {{ M[i, {c + 1}] = N[i, {k + 2}]; }}')
    if K:
        s = ' + '.join(f'N[i, {k + 2}]' for k in K)
        lines.append(f'for (i = 1; i <= 4; i++) {{ M[i, {len(cols) + 1}] = {s}; }}')
    eqK = ''.join(f', d{k} - sg' for k in K)
    lines += [f'ideal I = V{eqK}, minor(M, 4);',
              'poly opn = sg' + ''.join(f'*(d{k}^2 - sg^2)' for k in cols) + ';',
              'def S0 = sat(I, ideal(opn));',
              'ideal S; if (typeof(S0) == "list") { S = S0[1]; } else { S = S0; }',
              'S = std(S);',
              '"DIM"; dim(S);',
              'if (dim(S) > 0) { "DEG"; mult(S); }',
              'quit;']
    return '\n'.join(lines)


if __name__ == '__main__':
    reps, order = orbit_reps()
    print(f'tree group of order {order}; {len(reps)} orbits of zero patterns K')
    for K in reps:
        open('v5_tmp.sing', 'w').write(stratum_script(K))
        out = subprocess.run(['Singular', '-q', 'v5_tmp.sing'], capture_output=True, text=True, timeout=1800).stdout
        toks = [l.strip() for l in out.splitlines() if 'redefin' not in l and l.strip()]
        dim = toks[toks.index('DIM') + 1] if 'DIM' in toks else '?'
        deg = toks[toks.index('DEG') + 1] if 'DEG' in toks else '-'
        pdim = int(dim) - 1 if dim not in ('?',) and int(dim) > 0 else 'empty'
        print(f'  K = {K!s:<26} |K| = {len(K)}:  singular stratum projective dim {pdim}, degree {deg}', flush=True)
