"""Classify slice_search hits: nondegenerate E_4, identified against all known E_4."""
import sys
from collections import defaultdict
from slice import canonical, e4_at
from known_e4 import all_known, key, twin

known = {}
for name, r in all_known(kmax=14).items():
    known.setdefault(key(r), []).append(name)
    known.setdefault(key(twin(r)), []).append(name + '(twin)')
by = defaultdict(list)
deg = 0
for line in open(sys.argv[1]):
    P = tuple(int(v) for v in line.split())
    k = canonical(e4_at(*P))
    if k is None:
        deg += 1
        continue
    by[k].append(P)
print(f'{sum(len(v) for v in by.values())} nondegenerate points, {deg} degenerate, {len(by)} distinct E_4')
for k, pts in sorted(by.items(), key=lambda kv: kv[0][-1]):
    ids = known.get(k, ['NEW'])
    print(f'  max root {k[-1]:>40}  {ids[:3]}  {pts[:2]}')
