"""Degenerate slice points off the six split lines: fit low-degree curves."""
import sys
import itertools
import numpy as np
import sympy as sp
from slice import e4_at, canonical, x, y, z
from slice_search import search

LINES = [lambda a, b, c: c, lambda a, b, c: a + b, lambda a, b, c: a - c,
         lambda a, b, c: a + c, lambda a, b, c: b - c, lambda a, b, c: b + c]
hits = search(int(sys.argv[1]))
deg = [P for P in hits if canonical(e4_at(*P)) is None]
rest = [P for P in deg if all(f(*P) != 0 for f in LINES)]
print('degenerate points off the six lines:', len(rest))
# group by which |roots| coincide
from collections import Counter
pat = Counter()
for P in rest:
    r = [abs(v) for v in e4_at(*P)]
    eq = tuple((i, j) for i, j in itertools.combinations(range(8), 2) if r[i] == r[j])
    pat[eq] += 1
print('coincidence patterns:', pat.most_common(10))
# conic fit: monomials of degree 2; greedy grouping by exact nullspace of 5-point subsets
mons = lambda P: [P[0]**2, P[1]**2, P[2]**2, P[0]*P[1], P[0]*P[2], P[1]*P[2]]
pts = rest[:400]
found = {}
for sub in itertools.combinations(range(min(len(pts), 40)), 5):
    M = sp.Matrix([mons(pts[i]) for i in sub])
    ns = M.nullspace()
    if len(ns) != 1:
        continue
    v = ns[0]
    v = v / max(abs(c) for c in v if c != 0)
    den = sp.ilcm(*[sp.fraction(c)[1] for c in v])
    v = tuple(int(c * den) for c in v)
    g = 0
    for c in v:
        g = sp.igcd(g, c)
    v = tuple(c // g for c in v)
    if v in found:
        continue
    n = sum(1 for P in rest if sum(a * b for a, b in zip(v, mons(P))) == 0)
    if n >= 12:
        found[v] = n
for v, n in sorted(found.items(), key=lambda kv: -kv[1]):
    print('conic', sp.factor(sum(c * m for c, m in zip(v, [x**2, y**2, z**2, x*y, x*z, y*z]))), 'points', n)
