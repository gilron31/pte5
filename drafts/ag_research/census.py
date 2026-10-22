"""G3(f): which known E_4 coincide, and which small linear relations they satisfy.

A relation is c . r = 0 with c in {-1,0,1}^8 on the tree-ordered roots, tested on
every sign/order representative (signs are free; roots are only defined up to
sign) -- equivalently: is there a signed subset of the 8 |roots| summing to 0?
We record its support pattern per pair, e.g. (1,1,1,1) is Bremner's A-C=E-G
type and (2,2,2,2) the 8-term A+B+C+D=E+F+G+H type, and we test both twin forms.
"""
import itertools
from collections import defaultdict

from known_e4 import all_known, key, twin


def relations(r):
    """Support patterns (terms per pair, sorted within halves) of signed zero sums."""
    a = [abs(v) for v in r]
    found = set()
    for c in itertools.product((-1, 0, 1), repeat=8):
        nz = [i for i in range(8) if c[i]]
        if len(nz) < 3 or c[nz[0]] != 1:
            continue
        if sum(ci * ai for ci, ai in zip(c, a)) == 0:
            per_pair = tuple(sum(1 for i in (2*k, 2*k+1) if c[i]) for k in range(4))
            h1, h2 = sorted(per_pair[:2]), sorted(per_pair[2:])
            found.add(tuple(sorted([tuple(h1), tuple(h2)])))
    return found


if __name__ == '__main__':
    rows = all_known(kmax=12)
    ident = defaultdict(list)
    for name, r in rows.items():
        ident[key(r)].append((name, 'root form'))
        ident[key(twin(r))].append((name, 'twin form'))
    print('== coincidences between sources (same E_4 up to scaling) ==')
    for k, names in ident.items():
        srcs = {n for n, _ in names}
        if len(srcs) > 1:
            print('  ', names)
    print('\n== signed zero sums with >= 3 terms, as (terms per pair) patterns ==')
    for name, r in rows.items():
        rel = {('root',) + p for p in relations(r)} | {('twin',) + p for p in relations(twin(r))}
        print(f'  {name:>16}  {sorted(rel) if rel else "-"}')
