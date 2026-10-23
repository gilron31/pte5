"""N19: points of Sing(V_5) that N12/N12a missed, by Jacobian rank of V_5's 11
equations (7 pair-sum quadrics, then V_4's equations in d_k = r_{2k}^2 - r_{2k+1}^2).
Rank 11 = smooth; rank 7 = the corank of a cone-over-V_4 point (all four V_4-type
equations vanish to order >= 2).
  1. over the vertex of the base cone: all roots +-1 (exact, over Q);
  2. the "worst points" with other zero patterns, not only (1, 0)^8 (exact);
  3. sigma = 0 (r_{2k+1} = i r_{2k}): generic points, points with two pairs
     entirely zero (curves over a Fermat quartic), and points with one half of
     the tree zero (8 zero roots), over F_p.
Usage: python3 v5_audit_points.py
"""
import random

import sympy as sp

r = sp.symbols('r0:16')
d = [r[2*k]**2 - r[2*k+1]**2 for k in range(8)]
pair = [r[2*k]**2 + r[2*k+1]**2 for k in range(8)]
EQS = [pair[0] - pair[k] for k in range(1, 8)] + [
    d[0]**2 + d[1]**2 - d[2]**2 - d[3]**2, d[0]**2 + d[1]**2 - d[4]**2 - d[5]**2,
    d[0]**2 + d[1]**2 - d[6]**2 - d[7]**2, d[0]**2*d[1]**2 + d[2]**2*d[3]**2 - d[4]**2*d[5]**2 - d[6]**2*d[7]**2]
JAC = sp.lambdify(r, sp.Matrix(EQS).jacobian(r).tolist())
EV = sp.lambdify(r, EQS)
dd = sp.symbols('d0:8')
V4 = [dd[0]**2 + dd[1]**2 - dd[2]**2 - dd[3]**2, dd[0]**2 + dd[1]**2 - dd[4]**2 - dd[5]**2,
      dd[0]**2 + dd[1]**2 - dd[6]**2 - dd[7]**2, dd[0]**2*dd[1]**2 + dd[2]**2*dd[3]**2 - dd[4]**2*dd[5]**2 - dd[6]**2*dd[7]**2]
NJAC = sp.lambdify(dd, sp.Matrix(V4).jacobian(dd).tolist())


def rank_mod(M, p=None):
    if p is None:
        return sp.Matrix(M).rank()
    M = [[int(v) % p for v in row] for row in M]
    rk = 0
    for c in range(len(M[0])):
        piv = next((i for i in range(rk, len(M)) if M[i][c]), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][c], -1, p)
        M[rk] = [v * inv % p for v in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = M[i][c]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk


def rank_at(pt, p=None):
    vals = EV(*pt)
    assert all((v % p == 0) if p else v == 0 for v in vals), 'not on V5'
    return rank_mod(JAC(*pt), p)


def sigma0_point(x, i, p):
    out = []
    for v in x:
        out += [v % p, i * v % p]
    return out


if __name__ == '__main__':
    rng = random.Random(1)
    print('1. vertex points (all roots +-1):',
          [rank_at([rng.choice((1, -1)) for _ in range(16)]) for _ in range(4)])
    pts = []
    for _ in range(4):
        pt = []
        for k in range(8):
            a = rng.choice((1, -1))
            pt += [a, 0] if rng.random() < 0.5 else [0, a]
        pts.append(rank_at(pt))
    print('2. worst points, random zero patterns:', pts)
    # 3a. sigma = 0, generic: d = x^2 on V_4 over F_89
    p = 89
    i = next(a for a in range(2, p) if a * a % p == p - 1)
    fourth = {}
    for x in range(1, p):
        fourth.setdefault(pow(x, 4, p), []).append(x)
    gen = []
    while len(gen) < 6:
        x0, x1 = rng.randrange(1, p), rng.randrange(1, p)
        S = (pow(x0, 4, p) + pow(x1, 4, p)) % p
        sols = [(a, b) for a in range(1, p) for b in fourth.get((S - pow(a, 4, p)) % p, [])]
        if not sols:
            continue
        (x2, x3), (x4, x5) = rng.choice(sols), rng.choice(sols)
        T = (pow(x0 * x1, 4, p) + pow(x2 * x3, 4, p) - pow(x4 * x5, 4, p)) % p
        for a in range(1, p):
            b = (S - a) % p
            if a * b % p == T and a in fourth and b in fourth:
                x = [x0, x1, x2, x3, x4, x5, fourth[a][0], fourth[b][0]]
                dv = [v * v % p for v in x]
                gen.append((rank_at(sigma0_point(x, i, p), p), rank_mod(NJAC(*dv), p)))
                break
    print(f'3a. sigma = 0, all pairs nonzero, F_{p}: (V5 rank, V4 rank at d):', gen)
    # 3b. sigma = 0, pairs 0 and 1 zero: d3 = +-i d2, ..., d2^4 = d4^4 + d6^4 (Fermat quartic)
    sq_i = next(a for a in range(2, p) if a * a % p == i)
    two = []
    for x2 in range(1, p):
        for x4 in range(1, p):
            rhs = (pow(x2, 8, p) - pow(x4, 8, p)) % p
            x6 = next((v for v in range(1, p) if pow(v, 8, p) == rhs), None)
            if x6:
                two.append(rank_at(sigma0_point([0, 0, x2, sq_i * x2, x4, sq_i * x4, x6, sq_i * x6], i, p), p))
                break
        if len(two) >= 6:
            break
    print('3b. sigma = 0, two pairs zero:', two)
    # 3c. sigma = 0, pairs 0..3 zero: d = (0,0,0,0, 1, i, zeta8, i zeta8) up to scaling, over F_97
    p = 97
    g = next(a for a in range(2, p) if pow(a, 48, p) != 1)
    z16 = pow(g, 6, p)
    i, s_i = pow(z16, 4, p), pow(z16, 2, p)
    half = [rank_at(sigma0_point([0, 0, 0, 0, x4, s_i * x4, z16 * x4, s_i * z16 * x4], i, p), p) for x4 in range(1, 5)]
    print('3c. sigma = 0, a half of the tree zero:', half)
