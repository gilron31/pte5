"""Feed the best explicit Z_p-towers from local_upper.py to the tree-walk oracle."""
import sys
from fractions import Fraction

from local_upper import best_tower, seeds_from_catalog
from padic_oracle import Oracle, params_from_roots, check_tree

seeds = seeds_from_catalog(sys.argv[1] if len(sys.argv) > 1 else 'catalog_1e11.csv')
for p, ks in ((5, (5, 6)), (11, (5,)), (13, (3, 4)), (7, (11,))):
    cands = []
    for L1, form, s in seeds:
        b = best_tower(s, p)
        if b:
            cands.append((b, L1, form, s))
    (ht, top, vmu, L), L1, form, s = min(cands)
    mu = p ** vmu
    leaves = [x for r in s for x in (Fraction(L) + mu * r, Fraction(L) - mu * r)]
    ells = [Fraction(L)] + params_from_roots([mu * r for r in s])
    assert check_tree(s)
    for k in ks:
        if p ** k > 2_000_000:
            print(f'p={p} k={k}: skipped (table too large)')
            continue
        print(f'p={p} k={k}: claimed ht {ht}; oracle tuple floor {Oracle(p, k).floor(ells)} '
              f'(seed {L1} {form}, mu={mu}, L={L})')
