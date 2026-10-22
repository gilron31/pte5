"""G5: rank-1 Chabauty-Coleman for the genus-2 quotients certified by
g2_rank1_certify.gp, i.e. C_S: y^2 = prod_{e in S} (t + e) over Symes' E_4.

Goal: C_S(Q) = the 8 known points (the two points over t = oo, the trivial lift
mu = 0, and the 6 Weierstrass points t = -e, a zero root). If so, no E_5 lies
above Symes' E_4, for any scaling mu.

Model: the integral monic odd model W^2 = prod (X - A_i) (from g2descent.gp's
odd_model, then X = X'/D^2 clears denominators). Points: P0 = (0, W0) is the
image of t = oo; iota(P0) = (0, -W0); Weierstrass points (A_i, 0) and oo.

Steps, at a good prime p >= 5:
  1. I = int_{iota P0}^{P0} (dx/2y, x dx/2y) = 2 log(P0 - oo); omega = I1*w0 - I0*w1.
  2. For every F_p-point of the reduction, the residue disc D:
     - if D contains a known point Q (a Weierstrass point, oo, P0 or iota P0):
       int_Q omega over rational points vanishes (2-torsion or annihilated), so
       F(t) = tiny integral from Q in a local coordinate;
     - else F(t) = int_{P0}^{Q~} omega + tiny integral, with Q~ any lift.
     Count the zeros of F on the disc (t in pZ_p) by Strassmann / Newton polygon.
  3. Certificate iff each disc's zero count equals its number of known points.
Run (passagemath / Sage): python3 chabauty_g2.py
Tested with passagemath 10.8.13 (~/sage_venv); see chabauty_smoke.py for the A1 identity check.
"""
from sage.all import (QQ, ZZ, Qp, PolynomialRing, HyperellipticCurve, GF, lcm, prime_range,
                      PowerSeriesRing, Integer)

R_SYMES = [131, 343, 11, 367, 101, 353, 77, 359]
SIGNED = R_SYMES + [-r for r in R_SYMES]
CERTIFIED = [[1, 5, 6, 11, 12, 14], [6, 7, 8, 9, 14, 15], [1, 9, 11, 12, 14, 15],
             [2, 3, 8, 9, 11, 13], [3, 5, 8, 12, 13, 16], [5, 8, 9, 11, 14, 15], [3, 6, 7, 8, 9, 13]]


def odd_integral_model(e):
    e6 = QQ(e[5])
    c = 1
    for i in range(5):
        c *= (e6 - e[i])
    a = [-c / (e6 - e[i]) for i in range(5)]
    D = lcm([ai.denominator() for ai in a])
    A = [ZZ(ai * D**2) for ai in a]
    return A, D, c


def strassmann_zeros(series, p, N):
    """Upper bound for the zeros in pZ_p of F(t) = sum_k c_k t^k (power or Laurent
    series, indexed by exponent): with t = p s, the largest k minimizing v(c_k) + k
    (Strassmann). Only exponents k >= 0 occur for our integrals."""
    best_v, best_k = None, 0
    lo = max(series.valuation(), 0) if series != 0 else 0
    for k in range(0, min(N, series.prec() if series.prec() < 10**6 else N)):
        ck = series[k]
        if ck == 0:
            continue
        v = ck.valuation() + k
        if best_v is None or v < best_v:
            best_v, best_k = v, k
        elif v == best_v:
            best_k = k
    return best_k


def tiny_series(Cp, Q, omega_coeffs, prec):
    """Power series F(t) = int_Q^{Q(t)} (c0 + c1 x) dx/2y in a local coordinate t at Q."""
    x, y = Cp.local_coord(Q, prec)
    w = (omega_coeffs[0] + omega_coeffs[1] * x) * x.derivative() / (2 * y)
    return w.integral()


def chabauty(S, p=None, prec=20, tprec=40, verbose=True):
    e = [-SIGNED[i - 1] for i in S]
    A, D, c = odd_integral_model(e)
    Rx = PolynomialRing(QQ, 'x')
    xq = Rx.gen()
    f = 1
    for Ai in A:
        f *= (xq - Ai)
    disc = f.discriminant()
    W0sq = ZZ(f(0))
    if p is None:
        # good reduction, and the 8 known points in distinct residue discs (p does not divide f(0))
        p = next(q for q in prime_range(5, 1000) if ZZ(disc) % q != 0 and W0sq % q != 0)
    K = Qp(p, prec)
    C = HyperellipticCurve(f)
    Cp = C.change_ring(K)
    W0 = ZZ(f(0)).sqrt()
    P0, iP0 = Cp(0, W0), Cp(0, -W0)
    I = Cp.coleman_integrals_on_basis(iP0, P0)
    omega = (I[1], -I[0])                         # annihilates log J(Q)
    known_x = [0] + A                              # P0, iota P0 (x = 0), Weierstrass (A_i, 0), and oo
    Fp = GF(p)
    Cbar = HyperellipticCurve(f.change_ring(Fp))
    total_known = total_bound = 0
    report = []
    for Qbar in Cbar.points():
        if Qbar[2] == 0:                           # the point at infinity: disc of oo
            Q = Cp(0, 1, 0)
            known = 1
            F = tiny_series(Cp, Q, omega, tprec)
        else:
            xb, yb = ZZ(Qbar[0]), ZZ(Qbar[1])
            kn = [x0 for x0 in known_x if (x0 - xb) % p == 0]
            if yb == 0:                            # Weierstrass disc: base at the rational Weierstrass point
                x0 = next(a for a in A if (a - xb) % p == 0)
                Q = Cp(x0, 0)
                known = 1
                F = tiny_series(Cp, Q, omega, tprec)
            elif kn:                               # disc of P0 or iota P0
                Q = P0 if (W0 - yb) % p == 0 else iP0
                known = 1
                F = tiny_series(Cp, Q, omega, tprec)
            else:                                  # empty disc: lift a point, add the global Coleman integral
                Q = Cp.lift_x(K(xb))
                if (ZZ(Q[1].residue()) - yb) % p != 0:
                    Q = Cp(Q[0], -Q[1])
                J = Cp.coleman_integrals_on_basis(P0, Q)
                const = omega[0] * J[0] + omega[1] * J[1]
                F = const + tiny_series(Cp, Q, omega, tprec)
                known = 0
        bound = strassmann_zeros(F, p, tprec)
        total_known += known
        total_bound += bound
        report.append((tuple(Qbar), known, bound))
        if verbose and bound != known:
            print('   disc', tuple(Qbar), 'known', known, 'zero bound', bound)
    ok = total_bound == total_known
    if verbose:
        print(f'S = {S}, p = {p}: {len(report)} discs, known points {total_known}, '
              f'Chabauty bound {total_bound}  ->  {"C(Q) = known points: CERTIFIED" if ok else "not decided at this p"}')
    return ok, p, report


if __name__ == '__main__':
    for S in CERTIFIED:
        try:
            chabauty(S)
        except Exception as exc:  # report and continue with the next curve
            print(f'S = {S}: error {type(exc).__name__}: {exc}')
