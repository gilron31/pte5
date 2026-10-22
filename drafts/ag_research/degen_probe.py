"""Which curves carry the degenerate points of the slice search?"""
import sys
from collections import Counter
import sympy as sp
from slice import x, y, z, A_, B_, C_, D_, E_, G_, Q1, Q2, e4_at, canonical
from slice_search import search

hits = search(int(sys.argv[1]))
deg = [P for P in hits if canonical(e4_at(*P)) is None]
# test small candidate curves: lines and conics with small coefficients through many points
lines = Counter()
for (a, b, c) in deg:
    # record the reason: which roots vanish / coincide
    r = [abs(v) for v in e4_at(a, b, c)]
    zeros = tuple(i for i, v in enumerate(r) if v == 0)
    lines[('zeros', zeros)] += 1
print(lines.most_common(12))
# lines through the degenerate points: count incidences for lines ux+vy+wz=0 with |u|,|v|,|w|<=3
cnt = Counter()
rng = range(-3, 4)
for u in rng:
    for v in rng:
        for w in rng:
            if (u, v, w) <= (0, 0, 0) or (u, v, w) == (0, 0, 0):
                continue
            n = sum(1 for (a, b, c) in deg if u * a + v * b + w * c == 0)
            if n > 20:
                cnt[(u, v, w)] = n
print('lines with many degenerate points:', cnt.most_common(20))
print('total degenerate', len(deg))
