"""G5: Mordell-Weil sieve completing rank-1 Chabauty on a certified quotient C_S.

Facts used (rank J(Q) = 1, J(Q)_tors = J[2], D0 = [P0 - oo] of infinite order,
D0 not in 2J(Q) + J[2] -- all from g2_rank1_certify.gp):
 * Chabauty at p0 (chabauty_g2.py): every known disc contains exactly one p-adic
   zero, its known point; so an unknown rational point z lies in an EMPTY disc mod p0.
 * Saturation: for every odd prime l dividing #J(F_p) of a sieve prime p, D0 is
   shown not l-divisible mod torsion, using an auxiliary prime q with l || #J(F_q):
   there x in lJ(F_q) iff (#J(F_q)/l) x = 0. Then J(Q) = <G> + J[2] with D0 = mG + T0,
   gcd(m, ord_p) = 1, and the image of J(Q) in J(F_p) is <D0> + J[2].
 * So [z - oo] = n G + T reduces to c D0 + T mod p for a single integer
   c = n m^-1 (mod each ord_p) and a single T in J[2]. Sieve over (c mod L, T),
   L = lcm of ord_p(D0): require z mod p0 in an empty disc and z mod p on the curve.
If nothing survives, C(Q) = the 8 known points.
Usage: ~/sage_venv/bin/python mw_sieve.py [index of certified curve]
"""
import sys
from itertools import product

from sage.all import GF, ZZ, QQ, PolynomialRing, HyperellipticCurve, prime_range, lcm, divisors

from chabauty_g2 import odd_integral_model, SIGNED, CERTIFIED, chabauty


def setup(S):
    e = [-SIGNED[i - 1] for i in S]
    A, D, c = odd_integral_model(e)
    R = PolynomialRing(QQ, 'x')
    x = R.gen()
    f = 1
    for a in A:
        f *= (x - a)
    return f, A, ZZ(f(0)).sqrt()


class Jp:
    """J(F_p) of y^2 = f mod p with helpers."""
    def __init__(self, f, p):
        self.p = p
        self.C = HyperellipticCurve(f.change_ring(GF(p)))
        self.J = self.C.jacobian()
        self.X = self.J(GF(p))
        self.N = ZZ(self.C.frobenius_polynomial()(1))
        self.zero = self.X(0)

    def cls(self, P):
        """[P - oo] for an affine point P = (x, y) of the reduction."""
        return self.X(self.C(P[0], P[1]))

    def order(self, D):
        for d in divisors(self.N):
            if d * D == self.zero:
                return d
        raise ValueError('order not found')


def two_torsion(J, A):
    gens = [J.cls((a, 0)) for a in A]
    T = set()
    out = []
    for bits in product((0, 1), repeat=4):
        D = J.zero
        for b, g in zip(bits, gens[:4]):
            if b:
                D = D + g
        out.append(D)
    return out          # the 5th Weierstrass class is the sum of the other four


def saturated(f, A, W0, ell, bad, max_q=3000):
    """Prove D0 not ell-divisible mod J[2], via a prime q with ell || #J(F_q)."""
    for q in prime_range(5, max_q):
        if q in bad:
            continue
        Jq = Jp(f, q)
        if Jq.N % ell or (Jq.N // ell) % ell == 0:
            continue
        D0 = Jq.cls((0, W0))
        n1 = Jq.N // ell
        if all(n1 * (D0 - T) != Jq.zero for T in two_torsion(Jq, A)):
            return q
    return None


def sieve(idx, pool_max=120, max_steps=12, verbose=True):
    """Incremental MW sieve: survivors (c mod L, T) start from the Chabauty prime's
    empty-disc condition; primes are added greedily, choosing the one with the
    smallest expected growth (lift factor L'/L) * (#C(F_p) / #J(F_p))."""
    S = CERTIFIED[idx]
    f, A, W0 = setup(S)
    disc = ZZ(f.discriminant())
    bad = set(ZZ(disc * f(0)).prime_divisors()) | {2, 3}
    ok, p0, report = chabauty(S, verbose=False)
    assert all(b == k for (_, k, b) in report if k > 0), 'a known disc has extra p-adic zeros'
    empty = [Q for (Q, k, b) in report if k == 0]

    def prime_data(p, empty_pts=None):
        Jj = Jp(f, p)
        D0 = Jj.cls((0, W0))
        o = Jj.order(D0)
        tors = two_torsion(Jj, A)
        if empty_pts is None:
            target = {str(Jj.zero)}
            for Q in Jj.C.points():
                if Q[2] != 0:
                    target.add(str(Jj.cls((Q[0], Q[1]))))
        else:
            target = {str(Jj.cls((ZZ(Q[0]), ZZ(Q[1])))) for Q in empty_pts}
        mult = [(k * D0) for k in range(o)]
        ok_table = [[str(mult[k] + tors[t]) in target for t in range(16)] for k in range(o)]
        return dict(p=p, N=Jj.N, o=o, ok=ok_table, ncurve=len(target))

    d0 = prime_data(p0, empty)
    used = [d0]
    L = d0['o']
    surv = [(c, t) for c in range(L) for t in range(16) if d0['ok'][c][t]]
    if verbose:
        print(f'S = {S}; Chabauty prime p0 = {p0}: #J = {d0["N"]}, ord(D0) = {L}, '
              f'empty discs {len(empty)}; initial survivors {len(surv)}', flush=True)
    pool = {}
    for q in prime_range(5, pool_max):
        if q in bad or q == p0:
            continue
        pool[q] = prime_data(q)
    for step in range(max_steps):
        if not surv:
            break
        def growth(d):
            Ln = lcm(L, d['o'])
            return (Ln // L) * (d['ncurve'] / d['N'])
        q = min(pool, key=lambda q: growth(pool[q]))
        d = pool.pop(q)
        Ln = lcm(L, d['o'])
        new = []
        for c, t in surv:
            for k in range(Ln // L):
                cc = c + k * L
                if d['ok'][cc % d['o']][t]:
                    new.append((cc, t))
        surv, L = new, Ln
        used.append(d)
        if verbose:
            print(f'   + p = {q}: ord(D0) = {d["o"]}, modulus L = {L}, survivors {len(surv)}', flush=True)
    # saturation at the odd primes dividing the orders actually used
    ells = sorted({l for d in used for l in ZZ(d['o']).prime_divisors() if l != 2})
    sat_ok = True
    for l in ells:
        qq = saturated(f, A, W0, l, bad)
        sat_ok &= qq is not None
        if verbose:
            print(f'   saturation at l = {l}: ' + (f'proved via q = {qq}' if qq else 'NOT proved'), flush=True)
    verdict = (not surv) and sat_ok
    if verbose:
        print(f'S = {S}: survivors {len(surv)}, saturation {"ok" if sat_ok else "FAILED"} -> '
              + ('C(Q) = the 8 known points: CERTIFIED' if verdict else 'not certified'), flush=True)
    return verdict, surv


if __name__ == '__main__':
    idx = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    sieve(idx)
