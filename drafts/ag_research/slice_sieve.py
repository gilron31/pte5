"""R1 beyond int64: points of Bremner's slice with |x|, |y|, z <= H for large H.

For each small prime p, a table T_p[z mod p][x mod p][y mod p] marks residues
where q1 and q2 are both squares (or 0) mod p. A point survives only if it is
marked for every p; survivors (about 2^-2k of all points for k primes) are then
checked exactly with Python integers. Symmetries: z >= 1, x >= y.
Usage: python3 slice_sieve.py H [z_from z_to]
"""
import sys
from collections import defaultdict
from math import gcd, isqrt

import numpy as np

from slice import canonical, e4_at
from known_e4 import all_known, key, twin

PRIMES = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]


def q_mod(p, Z):
    r = np.arange(p, dtype=np.int64)
    X, Y = np.meshgrid(r, r, indexing='ij')
    out = []
    for z in (Z % p, (-Z) % p):
        X2, Y2 = X * X % p, Y * Y % p
        v = (X2 * X2 % p * Y2 - 2 * X2 * X % p * Y2 % p * Y + X2 * Y2 % p * Y2
             + z * z * (10 * X2 * X % p * Y - 20 * X2 * Y2 + 10 * X * Y2 % p * Y)
             + z**4 * (45 * X2 - 10 * X * Y + 45 * Y2)
             - 4 * z**3 * (X2 * X + X2 * Y + X * Y2 + Y2 * Y) + 40 * z**5 * (X + Y)) % p
        out.append(v)
    return out


def qr_table(p):
    sq = np.zeros(p, dtype=bool)
    sq[(np.arange(p) ** 2) % p] = True
    return sq


def search(H, z0, z1):
    rng = np.arange(-H, H + 1, dtype=np.int64)
    X, Y = np.meshgrid(rng, rng, indexing='ij')
    keep = (X >= Y) & (X + Y != 0)
    Xk, Yk = X[keep], Y[keep]
    res = {p: (Xk % p, Yk % p) for p in PRIMES}
    sq = {p: qr_table(p) for p in PRIMES}
    hits = []
    for Z in range(z0, z1 + 1):
        # drop the degenerate split lines x = +-z, y = +-z (and x + y = 0 above)
        mask = (np.abs(Xk) != Z) & (np.abs(Yk) != Z)
        for p in PRIMES:
            q1, q2 = q_mod(p, Z)
            ok = sq[p][q1] & sq[p][q2]
            xr, yr = res[p]
            mask &= ok[xr, yr]
        for i in np.nonzero(mask)[0]:
            a, b = int(Xk[i]), int(Yk[i])
            if gcd(gcd(a, b), Z) != 1:
                continue
            if e4_at(a, b, Z) is not None:
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
    known = {}
    for name, r in all_known(kmax=12).items():
        known.setdefault(key(r), []).append(name)
        known.setdefault(key(twin(r)), []).append(name + '(twin)')
    print(f'H={H} z in [{z0},{z1}]: {len(hits)} points, {len(by_e4)} nondegenerate E_4', flush=True)
    for k, pts in sorted(by_e4.items(), key=lambda kv: kv[0][-1]):
        print(f'  max root {k[-1]:>30}  identity {known.get(k, ["NEW"])}  points {pts[:4]}', flush=True)
