"""A5: run the bottom-lift certificate over every catalog E_4 (both 45-degree twin
forms). For each seed: the first random non-symmetric 6-subset whose quotient is
certified (rank 1, Chabauty, MW sieve) proves that no E_5 lies above that E_4,
for any scaling mu.
Usage: ~/sage_venv/bin/python -u certify_catalog.py [catalog.csv] [tries]
"""
import csv
import sys
import time
from math import isqrt

from sage.all import QQ

from lift_pipeline import certify_seed


def seeds(path):
    for row in csv.DictReader(open(path)):
        L1 = int(row['l1'])
        pts = []
        for i in range(4):
            nu = int(row[f'nu{i}'])
            pts.append((isqrt((L1 + nu) // 2), isqrt((L1 - nu) // 2)))
        yield L1, 'root', [v for m, n in pts for v in (m - n, m + n)]
        yield L1, 'point', [v for m, n in pts for v in (m, n)]


if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else 'catalog_1e12.csv'
    tries = int(sys.argv[2]) if len(sys.argv) > 2 else 80
    done = 0
    for L1, form, rho in seeds(path):
        t0 = time.time()
        res = certify_seed(rho, tries=tries, verbose=False)
        dt = time.time() - t0
        if res:
            done += 1
            S, p0 = res
            print(f'L1 = {L1:>13} {form:>5}: CERTIFIED no E_5 above (quotient {S}, p0 = {p0}; {dt:.0f} s)', flush=True)
        else:
            print(f'L1 = {L1:>13} {form:>5}: undecided after {tries} subsets ({dt:.0f} s)', flush=True)
    print(f'done: {done} seeds certified', flush=True)
