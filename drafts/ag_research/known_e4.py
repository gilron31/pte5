"""Every known E_4 in one place, as root vectors in tree order
(pairs (0,1),(2,3) | (4,5),(6,7)), each also in its 45-degree twin form.

Sources: the S-A catalog to 1e11 (drafts/ag_pilots/catalog_1e11.csv), Bremner
2008 Table 1 (rows 1-11, first line of each pairing), and the two families at
multiples k of the generator (code/old/utils.py, Bremner's parameterizations).
"""
import csv
import os
import sys
from fractions import Fraction
from math import gcd, isqrt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'code', 'old'))

# Bremner Table 1: A B C D | E F G H (pairs (A,B),(C,D),(E,F),(G,H)); * = on the slice
BREMNER_T1 = {
    '1*': (-101, 353, -77, 359, 343, 131, 367, 11),
    '2*': (-719, -827, 139, 1087, -541, 953, 317, -1049),
    '3*': (-141, 2428, 1267, 2076, 916, 2253, 2324, 717),
    '4': (5093, 9217, 9677, 4153, 8473, 6253, 3167, 10043),
    '5': (5971, 11287, 5447, 11549, 3661, 12233, 9689, 8317),
    '6': (4825, 17573, 8677, 16025, 9935, 15277, 17923, 3295),
    '7': (21523, 4567, 17593, 13213, 19483, 10223, 7493, 20687),
    '8': (3031, 32507, 30109, 12623, 7397, 31799, 10129, 31037),
    '9*': (-171061, 242103, -216799, 202173, -194223, 223949, -239961, 174053),
    '10*': (-50169, 2324657, 1368087, 1880129, 689137, 2220729, 2107393, 982569),
    '11*': (-42677861, 19377903, -46331907, 7089241, 15833813, 44115681, 12179767, 45260979),
}


def is_e4(r):
    a = [Fraction(v) ** 2 for v in r]
    n = {a[0] + a[1], a[2] + a[3], a[4] + a[5], a[6] + a[7]}
    return len(n) == 1 and a[0] * a[1] + a[2] * a[3] == a[4] * a[5] + a[6] * a[7]


def primitive(r):
    """Scale a rational root vector to coprime integers."""
    r = [Fraction(v) for v in r]
    den = 1
    for v in r:
        den = den * v.denominator // gcd(den, v.denominator)
    ints = [int(v * den) for v in r]
    g = 0
    for v in ints:
        g = gcd(g, v)
    return tuple(v // g for v in ints)


def twin(r):
    """The 45-degree twin: (A,B) -> (A+B, A-B) in every pair, made primitive."""
    return primitive([v for k in range(4) for v in (r[2*k] + r[2*k+1], r[2*k] - r[2*k+1])])


def key(r):
    """Shape-independent identity: sorted primitive |roots|."""
    return tuple(sorted(abs(v) for v in primitive(r)))


def catalog(path=os.path.join(HERE, '..', 'ag_pilots', 'catalog_1e11.csv')):
    out = {}
    for row in csv.DictReader(open(path)):
        L1 = int(row['l1'])
        pts = []
        for i in range(4):
            nu = int(row[f'nu{i}'])
            pts.append((isqrt((L1 + nu) // 2), isqrt((L1 - nu) // 2)))
        out[f'cat{L1}'] = primitive([v for m, n in pts for v in (m, n)])
    return out


def family(which, ks):
    import utils  # code/old/utils.py
    gen = utils.bremner_F1_point if which == 1 else utils.bremner_F2_point
    out = {}
    for k in ks:
        try:
            r = gen(k)
        except Exception:
            continue
        r = primitive(r)
        if 0 in r or len({abs(v) for v in r}) < 8:
            continue
        out[f'F{which}k{k}'] = r
    return out


def all_known(kmax=12):
    rows = {}
    rows.update({f'T1-{k}': primitive(v) for k, v in BREMNER_T1.items()})
    rows.update(catalog())
    rows.update(family(1, range(-kmax, kmax + 1)))
    rows.update(family(2, range(-kmax, kmax + 1)))
    for name, r in rows.items():
        assert is_e4(r), name
    return rows


if __name__ == '__main__':
    rows = all_known()
    print(len(rows), 'root vectors; all satisfy the E_4 equations')
