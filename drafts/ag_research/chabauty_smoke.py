"""A1 smoke test: Coleman integrals on a certified quotient must satisfy
int_W^{P0} omega = (1/2) int_{iota P0}^{P0} omega for each rational Weierstrass point W
(log(W - oo) = 0 since W - oo is 2-torsion)."""
from sage.all import QQ, ZZ, Qp, PolynomialRing, HyperellipticCurve, prime_range
from chabauty_g2 import odd_integral_model, SIGNED, CERTIFIED

S = CERTIFIED[0]
e = [-SIGNED[i - 1] for i in S]
A, D, c = odd_integral_model(e)
R = PolynomialRing(QQ, 'x'); x = R.gen()
f = 1
for a in A:
    f *= (x - a)
disc = ZZ(f.discriminant())
p = next(q for q in prime_range(5, 1000) if disc % q != 0)
print('S =', S, 'A =', A, 'D =', D, 'good prime p =', p)
K = Qp(p, 15)
C = HyperellipticCurve(f).change_ring(K)
W0 = ZZ(f(0)).sqrt()
P0, iP0 = C(0, W0), C(0, -W0)
I = C.coleman_integrals_on_basis(iP0, P0)
print('I = int_{iota P0}^{P0} (w0, w1, w2, w3) =', I)
for a in A:
    W = C(a, 0)
    try:
        J = C.coleman_integrals_on_basis(W, P0)
        print(f'  W = ({a}, 0): 2*J - I (holomorphic part) =', [2 * J[i] - I[i] for i in range(2)])
    except Exception as exc:
        print(f'  W = ({a}, 0): {type(exc).__name__}: {exc}')
