"""Local ceiling for top-node arguments: for each prime p, the smallest
v_p(sqrt L_5) over explicit Z_p-towers built from the catalog seeds.

A Z_p-tower with v_p(sqrt L_5) = t shows that no p-adic argument can force
p^(t+1) | sqrt(L_5). Compare with CMSV's exponent e_p = v_p(C'_16 / 2): if the
ceiling t equals e_p, CMSV is already locally optimal at p.
Usage: python3 top_node_witness.py [catalog.csv]
"""
import sys

sys.path.insert(0, '../ag_pilots')
from local_upper import seeds_from_catalog, seed_levels, vp, is_square_Zp  # noqa: E402

CMSV = {2: 10, 3: 6, 5: 4, 7: 3, 11: 2, 13: 2, 17: 1, 19: 1, 23: 1, 29: 1,
        31: 0, 37: 1, 41: 1, 43: 1, 47: 0, 53: 1}


def best_top(seeds, p, max_vmu=4, Lrange=200):
    best = None
    for L1, form, s in seeds:
        if p == 2 and form == 'point':
            continue            # half-integral nodes: keep p = 2 integral
        lev = seed_levels(s)
        for vmu in range(max_vmu + 1):
            mu = p ** vmu
            for L in range(1, Lrange):
                vals = [L + e * mu * r for r in s for e in (1, -1)]
                if 0 in vals or not all(is_square_Zp(a, p) for a in vals):
                    continue
                top = 8 * vmu + vp(lev[3][0], p)
                ht = max([max(vp(a, p) for a in vals) // 2] +
                         [2 ** (j - 1) * vmu + max(vp(t, p) for t in lev[j - 1]) for j in range(1, 5)])
                cand = (top, ht, L1, form, vmu, L)
                if best is None or cand < best:
                    best = cand
                break                       # smallest L for this v(mu) suffices
            if best is not None and best[0] <= 8 * vmu:
                break
    return best


if __name__ == '__main__':
    seeds = seeds_from_catalog(sys.argv[1] if len(sys.argv) > 1 else 'catalog_1e12.csv')
    print(' p   CMSV e_p   witness top-node ceiling   verdict              [ht, seed, v(mu), L]')
    for p, e in CMSV.items():
        top, ht, L1, form, vmu, L = best_top(seeds, p)
        verdict = 'CMSV locally optimal' if top <= e else f'room: up to p^{top}'
        print(f'{p:>3}   {e:>4}       {top:>4}                     {verdict:<20} [{ht}, {L1} {form}, {vmu}, {L}]')
