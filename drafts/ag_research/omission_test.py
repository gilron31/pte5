"""Soundness test of the Chabauty + MW-sieve certificate by omission: hide a genuine
rational point (treat its disc mod p0 as empty); the sieve must keep its class
(c mod L, T) among the survivors and must not certify."""
from sage.all import QQ, ZZ, GF
from lift_pipeline import odd_model, poly_from, chabauty, mw_sieve, Jp, two_torsion

R = [131, 343, 11, 367, 101, 353, 77, 359]
signed = R + [-r for r in R]
S = [1, 5, 6, 11, 12, 14]                        # certified curve 0
e = [-signed[i - 1] for i in S]
A, to_X = odd_model(e)
f = poly_from(A)
W0 = QQ(f(0)).sqrt()
ok, p0, report = chabauty(A, [QQ(0)])
print('Chabauty prime', p0, 'known discs fine:', ok)
J0 = Jp(f, p0)
D0 = J0.cls(0, W0)
tors = two_torsion(J0, A)
cases = {
    'iota P0 (class -D0)': ((0, -W0), None),
    'Weierstrass (A1, 0)': ((A[0], 0), None),
    'Weierstrass (A3, 0)': ((A[2], 0), None),
}
for name, ((x, y), _) in cases.items():
    hide = [(GF(p0)(x), GF(p0)(y), 1)]
    cert, surv, log, (L, used) = mw_sieve(A, [QQ(0)], p0, report, hide=hide)
    # the true class of the hidden point, as (c, T): find it by brute force mod the first prime
    target = J0.cls(x, y)
    true_ct = [(c, t) for c in range(used[0]['o']) for t in range(16) if c * D0 + tors[t] == target]
    kept = [st for st in surv if (st[0] % used[0]['o'], st[1]) in true_ct]
    print(f'hide {name}: certified = {cert}; survivors {len(surv)}; hidden class kept: {bool(kept)}')
    print('    ' + '; '.join(log))
