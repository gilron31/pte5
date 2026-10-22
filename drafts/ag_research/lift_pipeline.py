"""G5 pipeline: certify that no lift exists above a seed, for any scaling mu.

Seed: positive roots rho of an E_{d-1} (4 roots for an E_3 seed, 8 for an E_4).
A lift is a t = L/mu with every t +- rho_k in one square class. Any 6 of the 2m
signed values give a genus-2 quotient C_e: y^2 = prod_{i} (t - e_i) containing
the lift as a rational point. Known rational points of every C_e: t = oo (mu = 0)
and t = e_i (a zero root) -- all degenerate -- plus any extra known lifts.

Per quotient:
  1. rank J(Q) = 1 certificate (g2_rank1_generic.gp: 2-descent, torsion, infinite order);
  2. Chabauty at a good prime p0: every disc holding a known point must contain
     exactly that one p-adic zero (Strassmann);
  3. Mordell-Weil sieve: no (c mod L, T) with the class c*D0 + T landing in an
     empty disc mod p0 and on the curve mod the sieve primes; saturation of D0 at
     the odd primes dividing the orders used.
Success on one quotient: the rational points of C_e are exactly the known ones,
so no lift beyond the known ones exists above the seed.
Usage (Sage): see the __main__ block.
"""
import random
import subprocess
from itertools import product

from sage.all import (QQ, ZZ, Qp, GF, PolynomialRing, HyperellipticCurve, prime_range, lcm,
                      divisors)


# ---------------------------------------------------------------- models
def odd_model(e):
    """t = e6 + 1/s, X = c s; then X' = D^2 X: W'^2 = prod (X' - A_i), integral monic.
    Returns A (5 integers), and the map t -> X'."""
    e = [QQ(v) for v in e]
    e6 = e[5]
    c = 1
    for i in range(5):
        c *= (e6 - e[i])
    a = [-c / (e6 - e[i]) for i in range(5)]
    D = lcm([ai.denominator() for ai in a])
    A = [ZZ(ai * D**2) for ai in a]
    to_X = lambda t: None if t is None else (QQ(0) if t == 'oo' else D**2 * c / (QQ(t) - e6))
    return A, to_X


def poly_from(A):
    R = PolynomialRing(QQ, 'x')
    x = R.gen()
    f = 1
    for a in A:
        f *= (x - a)
    return f


# ---------------------------------------------------------------- step 1: rank 1
def rank1(e):
    cmd = ('default(parisizemax, 2^31);\nread("g2_rank1_generic.gp");\n'
           f'print(certify_e({[str(QQ(v)) for v in e]}));\n').replace("'", '')
    out = subprocess.run(['gp', '-q'], input=cmd, capture_output=True, text=True, timeout=600).stdout
    line = [l for l in out.splitlines() if l.startswith('[')]
    return bool(line) and line[-1].startswith('[1')


# ---------------------------------------------------------------- step 2: Chabauty
def strassmann(F, N):
    best_v, best_k = None, 0
    top = N
    try:
        top = min(N, int(F.prec()))
    except (TypeError, ValueError, OverflowError):
        pass
    for k in range(0, top):
        ck = F[k]
        if ck == 0:
            continue
        v = ck.valuation() + k
        if best_v is None or v < best_v:
            best_v, best_k = v, k
        elif v == best_v:
            best_k = k
    return best_k


def tiny(Cp, Q, om, prec):
    x, y = Cp.local_coord(Q, prec)
    return ((om[0] + om[1] * x) * x.derivative() / (2 * y)).integral()


def chabauty(A, known_X, p=None, prec=20, tprec=40):
    """known_X: X'-coordinates of known non-Weierstrass rational points (each with both
    signs of W'), including 0 for t = oo. Returns (ok, p, report, Cp, P0)."""
    f = poly_from(A)
    disc = ZZ(f.discriminant())
    knownW = []
    for X in known_X:
        w2 = QQ(f(X))
        assert w2.is_square(), 'known point not on the curve'
        knownW.append((QQ(X), w2.sqrt()))
    def separated(q):
        if disc % q == 0:
            return False
        if any(QQ(X).denominator() % q == 0 for X, _ in knownW):
            return False
        # known points (incl. Weierstrass) in distinct residue discs
        res = []
        for X, W in knownW:
            res += [(GF(q)(X), GF(q)(W)), (GF(q)(X), GF(q)(-W))]
        res += [(GF(q)(a), GF(q)(0)) for a in A]
        return len(set(res)) == len(res)
    if p is None:
        p = next(q for q in prime_range(5, 2000) if separated(q))
    K = Qp(p, prec)
    Cp = HyperellipticCurve(f).change_ring(K)
    X0, W0 = knownW[0]                                   # t = oo: X' = 0
    P0, iP0 = Cp(K(X0), K(W0)), Cp(K(X0), -K(W0))
    I = Cp.coleman_integrals_on_basis(iP0, P0)
    om = (I[1], -I[0])
    Cbar = HyperellipticCurve(f.change_ring(GF(p)))
    report = []
    for Qb in Cbar.points():
        if Qb[2] == 0:
            Q, known = Cp(0, 1, 0), 1
            F = tiny(Cp, Q, om, tprec)
        else:
            xb, yb = ZZ(Qb[0]), ZZ(Qb[1])
            if yb == 0:
                a = next(a for a in A if (a - xb) % p == 0)
                Q, known = Cp(a, 0), 1
                F = tiny(Cp, Q, om, tprec)
            else:
                hit = [(X, s * W) for X, W in knownW for s in (1, -1)
                       if GF(p)(X) == GF(p)(xb) and GF(p)(s * W) == GF(p)(yb)]
                if hit:
                    X, W = hit[0]
                    Q, known = Cp(K(X), K(W)), 1
                    F = tiny(Cp, Q, om, tprec)
                else:
                    Q = Cp.lift_x(K(xb))
                    if GF(p)(ZZ(Q[1].residue())) != GF(p)(yb):
                        Q = Cp(Q[0], -Q[1])
                    J = Cp.coleman_integrals_on_basis(P0, Q)
                    F = om[0] * J[0] + om[1] * J[1] + tiny(Cp, Q, om, tprec)
                    known = 0
        report.append((tuple(Qb), known, strassmann(F, tprec)))
    ok = all(b == k for (_, k, b) in report if k > 0)
    return ok, p, report


# ---------------------------------------------------------------- step 3: MW sieve
class Jp:
    def __init__(self, f, p):
        self.p = p
        self.C = HyperellipticCurve(f.change_ring(GF(p)))
        self.X = self.C.jacobian()(GF(p))
        self.N = ZZ(self.C.frobenius_polynomial()(1))
        self.zero = self.X(0)

    def cls(self, x, y):
        return self.X(self.C(x, y))

    def order(self, D):
        for d in divisors(self.N):
            if d * D == self.zero:
                return d


def two_torsion(Jj, A):
    gens = [Jj.cls(a, 0) for a in A[:4]]
    out = []
    for bits in product((0, 1), repeat=4):
        D = Jj.zero
        for b, g in zip(bits, gens):
            if b:
                D = D + g
        out.append(D)
    return out


def saturated(f, A, W0, ell, bad, max_q=3000):
    for q in prime_range(5, max_q):
        if q in bad:
            continue
        Jq = Jp(f, q)
        if Jq.N % ell:
            continue
        D0 = Jq.cls(0, W0)
        n1 = Jq.N // ell
        if all(n1 * (D0 - T) != Jq.zero for T in two_torsion(Jq, A)):
            return q
    return None


def mw_sieve(A, known_X, p0, report, pool_max=120, max_steps=14, hide=()):
    f = poly_from(A)
    W0 = QQ(f(0)).sqrt()
    disc = ZZ(f.discriminant())
    bad = set(ZZ(disc * f(0)).prime_divisors()) | {2, 3}
    for X in known_X:
        bad |= set(ZZ(QQ(X).denominator()).prime_divisors())
    empty = [Q for (Q, k, b) in report if k == 0] + list(hide)

    def data(p, empty_pts=None):
        Jj = Jp(f, p)
        D0 = Jj.cls(0, W0)
        o = Jj.order(D0)
        tors = two_torsion(Jj, A)
        if empty_pts is None:
            target = {str(Jj.zero)} | {str(Jj.cls(Q[0], Q[1])) for Q in Jj.C.points() if Q[2] != 0}
        else:
            target = {str(Jj.cls(ZZ(Q[0]), ZZ(Q[1]))) for Q in empty_pts}
        mult = [k * D0 for k in range(o)]
        return dict(p=p, N=Jj.N, o=o, ncurve=len(target),
                    ok=[[str(mult[k] + tors[t]) in target for t in range(16)] for k in range(o)])

    d0 = data(p0, empty)
    L, used = d0['o'], [d0]
    surv = [(c, t) for c in range(L) for t in range(16) if d0['ok'][c][t]]
    log = [f'p0 = {p0}: ord {L}, empty discs {len(empty)}, survivors {len(surv)}']
    pool = {q: data(q) for q in prime_range(5, pool_max) if q not in bad and q != p0}
    for _ in range(max_steps):
        if not surv or not pool or (hide and len(surv) <= 3 * len(hide)):
            break
        q = min(pool, key=lambda q: (lcm(L, pool[q]['o']) // L) * pool[q]['ncurve'] / pool[q]['N'])
        d = pool.pop(q)
        Ln = lcm(L, d['o'])
        surv = [(c + k * L, t) for c, t in surv for k in range(Ln // L)
                if d['ok'][(c + k * L) % d['o']][t]]
        L = Ln
        used.append(d)
        log.append(f'+ {q}: ord {d["o"]}, L = {L}, survivors {len(surv)}')
    ells = sorted({l for d in used for l in ZZ(d['o']).prime_divisors() if l != 2})
    sat = {l: saturated(f, A, W0, l, bad) for l in ells}
    log.append('saturation: ' + ', '.join(f'{l} via {q}' if q else f'{l} FAILED' for l, q in sat.items()))
    return (not surv) and all(sat.values()), surv, log, (L, used)


# ---------------------------------------------------------------- driver
def certify_seed(rho, extra_t=(), tries=60, seed=1, verbose=True):
    """Try random non-symmetric 6-subsets of the signed roots until one quotient is
    certified. extra_t: known lift values t (rational), assumed to be lifts."""
    signed = [QQ(r) for r in rho] + [-QQ(r) for r in rho]
    rng = random.Random(seed)
    for attempt in range(tries):
        S = sorted(rng.sample(range(len(signed)), 6))
        vals = [signed[i] for i in S]
        if set(vals) == set(-v for v in vals):
            continue
        e = [-v for v in vals]                         # C_e: y^2 = prod (t + signed) = prod (t - e)
        if not rank1(e):
            continue
        A, to_X = odd_model(e)
        known_X = [QQ(0)] + [to_X(t) for t in extra_t if QQ(t) != e[5]]
        ok, p0, report = chabauty(A, known_X)
        if not ok:
            if verbose:
                print(f'  S = {S}: Chabauty at {p0} leaves an extra zero in a known disc; next')
            continue
        cert, surv, log, _ = mw_sieve(A, known_X, p0, report)
        if verbose:
            print(f'  S = {S} (roots {[str(v) for v in vals]}): ' + '; '.join(log))
        if cert:
            return S, p0
        if verbose:
            print(f'     NOT certified: {len(surv)} survivors')
    return None


if __name__ == '__main__':
    import sys
    which = sys.argv[1] if len(sys.argv) > 1 else 'calib'
    if which == 'calib':
        # A3: Symes' derived E_3, with the known lift t0 = 67405/12 (and -t0)
        rho = [4187, 5607, 4767, 5123]
        print('calibration WITHOUT the known lift (must NOT certify):')
        print(' ->', certify_seed(rho, extra_t=(), tries=20))
        print('calibration WITH t = +-67405/12 known (should certify unless more lifts exist):')
        print(' ->', certify_seed(rho, extra_t=(QQ(67405) / 12, -QQ(67405) / 12), tries=20))
