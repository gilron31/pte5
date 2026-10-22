"""Multi-modular implicitization of a family's (J1, J2)-image curve.

The elimination over Q (implicit_fam2_stu.sing) stalls on coefficient growth,
but over F_p it takes about a second. So: eliminate over F_p for primes just
below 2^31, normalize the generator to leading coefficient 1, CRT the
coefficients, rationally reconstruct, and stop when two consecutive
reconstructions agree. Then verify at fresh primes that the result vanishes on
the curve (F_p-points of the (s:t:u)-curve map into the reconstructed curve).
Usage: python3 modular_image.py fam2
"""
import re
import subprocess
import sys
from fractions import Fraction

import sympy as sp

from phi_degree import invariants
from twin_check import CURVES

s, t, u, J1s, J2s = sp.symbols('s t u J1 J2')


def primes_below(n, k):
    out, p = [], n
    while len(out) < k:
        p -= 1
        if sp.isprime(p):
            out.append(p)
    return out


def script(curve, p):
    L1, L2, L3, _ = invariants(1, t, u)
    G = sp.expand(sp.sympify(CURVES[curve]).subs(s, 1))
    f = lambda e: str(sp.expand(e)).replace('**', '^')
    return f'''ring r = {p}, (T,t,u,J1,J2), dp;
short = 0;
poly L1 = {f(L1)};
poly L2 = {f(L2)};
poly L3 = {f(L3)};
poly G = {f(G)};
ideal I = G, T*L1 - 1, J1*L1^2 - L2, J2*L1^4 - L3;
ideal K = eliminate(I, T*t*u);
if (size(K) != 1) {{ "BAD"; size(K); quit; }}
poly F = K[1] / leadcoef(K[1]);
print(F);
quit;
'''


def parse(text):
    """Polynomial in J1, J2 printed with short = 0 -> {(a, b): coeff mod p}."""
    terms = {}
    expr = text.strip().replace('\n', '').replace(' ', '')
    for m in re.finditer(r'([+-]?)([0-9]*)\*?((?:J[12](?:\^[0-9]+)?\*?)*)', expr):
        sign, coef, mon = m.groups()
        if not coef and not mon:
            continue
        c = int(coef) if coef else 1
        if sign == '-':
            c = -c
        a = b = 0
        for v, e in re.findall(r'J([12])(?:\^([0-9]+))?', mon):
            if v == '1':
                a = int(e) if e else 1
            else:
                b = int(e) if e else 1
        terms[(a, b)] = terms.get((a, b), 0) + c
    return terms


def ratrec(a, m):
    """Rational reconstruction of a mod m (Wang), or None."""
    a %= m
    bound = int((m // 2) ** 0.5)
    r0, r1, s0, s1 = m, a, 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)


def reconstruct(curve, max_primes=40):
    primes = primes_below(2**31 - 1, max_primes + 4)
    acc, M, prev = None, 1, None
    for i, p in enumerate(primes[:max_primes]):
        open('modimg_tmp.sing', 'w').write(script(curve, p))
        out = subprocess.run(['Singular', '-q', 'modimg_tmp.sing'], capture_output=True, text=True, timeout=1200).stdout
        out = '\n'.join(l for l in out.splitlines() if 'redefin' not in l)
        if 'BAD' in out:
            print('prime', p, 'bad: elimination ideal not principal'); continue
        terms = {k: v % p for k, v in parse(out).items()}
        if acc is None:
            acc, M = dict(terms), p
        else:
            if set(terms) != set(acc):
                print('prime', p, 'unlucky (support differs), skipped'); continue
            for k in acc:      # CRT
                a, b = acc[k], terms[k]
                acc[k] = (a + M * ((b - a) * pow(M, -1, p) % p)) % (M * p)
            M *= p
        rec = {k: ratrec(v, M) for k, v in acc.items()}
        ok = all(v is not None for v in rec.values())
        print(f'prime {i + 1}: modulus {M.bit_length()} bits, {len(acc)} terms, reconstruction {"ok" if ok else "not yet"}', flush=True)
        if ok and rec == prev:
            return rec, primes[max_primes:]
        prev = rec if ok else None
    return None, []


def verify(curve, F, fresh):
    """At fresh primes p: F(J(P)) = 0 for every F_p-point P of the curve with L1(P) != 0."""
    L1, L2, L3, _ = invariants(1, t, u)
    G = sp.Poly(sp.expand(sp.sympify(CURVES[curve]).subs(s, 1)), t, u)
    fl = [sp.lambdify((t, u), e, 'math') for e in (L1, L2, L3)]
    for p in fresh[:2]:
        p = int(p)
        Fp = {k: v.numerator * pow(v.denominator, -1, p) % p for k, v in F.items()}
        checked = bad = 0
        for T_ in range(0, 400):
            g = sp.Poly(G.as_expr().subs(t, T_), u, modulus=p)
            for U_ in range(0, 400):
                if g.eval(U_) % p:
                    continue
                a, b2, c = (int(f(T_, U_)) % p for f in fl)
                if a == 0:
                    continue
                j1 = b2 * pow(a * a, -1, p) % p
                j2 = c * pow(pow(a, 4, p), -1, p) % p
                val = sum(cf * pow(j1, i, p) * pow(j2, j, p) for (i, j), cf in Fp.items()) % p
                checked += 1
                bad += val != 0
        print(f'verify at p = {p}: {checked} curve points checked, {bad} failures')


if __name__ == '__main__':
    curve = sys.argv[1]
    F, fresh = reconstruct(curve)
    if F is None:
        print('no stable reconstruction'); sys.exit(1)
    den = sp.ilcm(*[v.denominator for v in F.values()])
    ints = {k: int(v * den) for k, v in F.items()}
    g = 0
    for v in ints.values():
        g = sp.igcd(g, v)
    ints = {k: v // g for k, v in ints.items()}
    poly = ' + '.join(f'({v})*J1^{a}*J2^{b}' for (a, b), v in sorted(ints.items()))
    open(f'{curve}_F_mod.txt', 'w').write(poly)
    print(f'{curve} image: {len(ints)} terms, max coefficient {max(abs(v) for v in ints.values()).bit_length()} bits; '
          f'degrees J1 {max(a for a, _ in ints)}, J2 {max(b for _, b in ints)}')
    verify(curve, F, fresh)
