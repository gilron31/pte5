"""R1: rational points of small height on Bremner's slice, classified by curve.

Vectorized over (x, y) for each z; q1, q2 stay below 2^53 for H <= 170, so the
float square root followed by an exact int64 check is exact.
Usage: python3 slice_search.py H
"""
import sys
from collections import defaultdict
from math import gcd

import numpy as np

from slice import canonical, e4_at, on_curves


def q_np(X, Y, Z, sgn):
    z = sgn * Z
    return (X**4 * Y**2 - 2 * X**3 * Y**3 + X**2 * Y**4 + 10 * X**3 * Y * z**2
            - 20 * X**2 * Y**2 * z**2 + 10 * X * Y**3 * z**2 - 4 * X**3 * z**3
            - 4 * X**2 * Y * z**3 - 4 * X * Y**2 * z**3 - 4 * Y**3 * z**3
            + 45 * X**2 * z**4 - 10 * X * Y * z**4 + 45 * Y**2 * z**4
            + 40 * X * z**5 + 40 * Y * z**5)


def is_sq(v):
    r = np.round(np.sqrt(np.maximum(v, 0).astype(np.float64))).astype(np.int64)
    return (v >= 0) & (r * r == v)


def search(Hmax):
    rng = np.arange(-Hmax, Hmax + 1, dtype=np.int64)
    X, Y = np.meshgrid(rng, rng, indexing='ij')
    hits = []
    for Z in range(1, Hmax + 1):
        ok = is_sq(q_np(X, Y, Z, 1)) & is_sq(q_np(X, Y, Z, -1))
        for i, j in zip(*np.nonzero(ok)):
            a, b = int(X[i, j]), int(Y[i, j])
            if gcd(gcd(a, b), Z) == 1:
                hits.append((a, b, Z))
    return hits


if __name__ == '__main__':
    Hmax = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    hits = search(Hmax)
    by_e4 = defaultdict(list)
    degenerate = 0
    for P in hits:
        key = canonical(e4_at(*P))
        if key is None:
            degenerate += 1
            continue
        by_e4[key].append(P)
    print(f'H = {Hmax}: {len(hits)} primitive points with q1, q2 squares; '
          f'{degenerate} degenerate; {len(by_e4)} distinct nondegenerate E_4')
    tally = defaultdict(int)
    for key, pts in sorted(by_e4.items(), key=lambda kv: kv[0][-1]):
        curves = sorted({c for P in pts for c in on_curves(*P)})
        tally[tuple(curves)] += 1
        print(f'  max root {key[-1]:>22}  curves {curves or ["NONE"]}  points {pts[:3]}')
    print('tally by curve set:', dict(tally))
