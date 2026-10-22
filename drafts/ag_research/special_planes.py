"""R4: all planes G = aA + bC + cE in the ACEG model of V_4 on which the
denominator den = A^2 + C^2 - E^2 - G^2 divides the numerator N of
b = B^2 = N / den, so that B^2, D^2, F^2, H^2 all become quadrics there
(Bremner's slice A - C = E - G is a = -1, b = 1, c = 1).

Method: substitute G, divide N by den as polynomials in E (den is quadratic in
E when c^2 != -1), and demand that the remainder vanish identically in (A, C, E).
"""
import sympy as sp

A, C, E, a, b, c = sp.symbols('A C E a b c')
G = a*A + b*C + c*E
den = sp.expand(A**2 + C**2 - E**2 - G**2)
N = sp.expand(A**2*(E**2 + G**2 - C**2) + C**4 - E**4 - G**4)

q, r = sp.div(sp.Poly(N, E), sp.Poly(den, E))
eqs = [sp.expand(co) for co in sp.Poly(r.as_expr(), A, C, E).coeffs()]
eqs = [sp.numer(sp.together(e)) for e in eqs]
print('conditions:', len(eqs))
sols = sp.solve(eqs, [a, b, c], dict=True)
real = []
for s in sols:
    if all(v.is_real for v in s.values()) and len(s) == 3:
        real.append(s)
    print(' ', s)
print('real solutions with all of a, b, c determined:', len(real))
for s in real:
    Gs = G.subs(s)
    print('   G =', Gs, '  den on plane =', sp.factor(den.subs(s)),
          '  B^2 on plane =', sp.factor(sp.cancel(N.subs(s) / den.subs(s))))
