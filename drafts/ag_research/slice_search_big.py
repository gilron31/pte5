"""R1 at larger height: exact int64 search of Bremner's slice up to H <= 500.

|q| <= 240 H^6 < 2^63 for H <= 500, and every intermediate term is smaller, so
the polynomial is evaluated exactly in int64.  Squares are tested by a float
estimate corrected by +-2 with exact int64 squares.  Symmetries used: z >= 1
(projective sign) and x >= y (q1, q2 are symmetric in x, y).
Usage: python3 slice_search_big.py H [z_from z_to]
"""
import sys
from collections import defaultdict
from math import gcd

import numpy as np

from slice import canonical, e4_at, on_curves
from known_e4 import all_known, key, twin


def q_pair(X, Y, Z):
    """(q1, q2) at the grid (X, Y) and integer Z, exactly in int64."""
    X2, Y2 = X * X, Y * Y
    even = (X2 * X2 * Y2 - 2 * X2 * X * Y2 * Y + X2 * Y2 * Y2
            + Z * Z * (10 * X2 * X * Y - 20 * X2 * Y2 + 10 * X * Y2 * Y)
            + Z**4 * (45 * X2 - 10 * X * Y + 45 * Y2))
    odd = (-4 * Z**3 * (X2 * X + X2 * Y + X * Y2 + Y2 * Y) + 40 * Z**5 * (X + Y))
    return even + odd, even - odd


def is_sq(v):
    ok = v >= 0
    r = np.floor(np.sqrt(np.where(ok, v, 0).astype(np.float64))).astype(np.int64)
    hit = np.zeros(v.shape, dtype=bool)
    for d in (-2, -1, 0, 1, 2):
        rr = r + d
        hit |= (rr >= 0) & (rr * rr == v)
    return ok & hit


def search(H, z0=1, z1=None):
    z1 = H if z1 is None else z1
    rng = np.arange(-H, H + 1, dtype=np.int64)
    X, Y = np.meshgrid(rng, rng, indexing='ij')
    keep = X >= Y
    Xk, Yk = X[keep], Y[keep]
    hits = []
    for Z in range(z0, z1 + 1):
        q1, q2 = q_pair(Xk, Yk, np.int64(Z))
        ok = is_sq(q1)
        idx = np.nonzero(ok)[0]
        if idx.size:
            ok2 = is_sq(q2[idx])
            for i in idx[ok2]:
                a, b = int(Xk[i]), int(Yk[i])
                if gcd(gcd(a, b), Z) == 1:
                    hits.append((a, b, Z))
    return hits


if __name__ == '__main__':
    H = int(sys.argv[1])
    z0, z1 = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (1, H)
    hits = search(H, z0, z1)
    by_e4 = defaultdict(list)
    for P in hits:
        k = canonical(e4_at(*P))
        if k is not None:
            by_e4[k].append(P)
    print(f'H={H} z in [{z0},{z1}]: {len(hits)} points, {len(by_e4)} nondegenerate E_4')
    known = {}
    for name, r in all_known(kmax=12).items():
        known.setdefault(key(r), []).append(name)
        known.setdefault(key(twin(r)), []).append(name + '(twin)')
    for k, pts in sorted(by_e4.items(), key=lambda kv: kv[0][-1]):
        names = known.get(k, ['NEW'])
        print(f'  max root {k[-1]:>26}  identity {names}  points {pts[:4]}')
