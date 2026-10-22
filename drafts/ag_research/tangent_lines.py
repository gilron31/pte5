"""Lines in Bremner's (s:t:u)-plane tangent to several of the branch conics.

Over a line that is tangent to k_t of the conics b, d, f, h (and meets the rest
transversally in distinct points), the (Z/2)^4-cover restricts to a cover with
4 - k_t nontrivial square conditions: genus 0 if k_t = 3, genus 1 if k_t = 2.
Such a line defined over Q, with the tangent restrictions *rational* squares
(not c * square with c a non-square), carries a genus-0/1 curve of S_B.
"""
import itertools

import sympy as sp

s, t, u, a, b_, c = sp.symbols('s t u a b c')
CON = {'b': 5*s**2 + 2*s*t + u**2, 'd': 5*s**2 - 2*s*t + u**2,
       'f': 5*s**2 + 2*s*u + t**2, 'h': 5*s**2 - 2*s*u + t**2}


def matrix(q):
    return sp.Matrix(3, 3, lambda i, j: sp.Rational(1, 2) * sp.diff(q, [s, t, u][i], [s, t, u][j]))


def dual(q):
    l = sp.Matrix([a, b_, c])
    return sp.expand((l.T * matrix(q).adjugate() * l)[0])


DUAL = {k: dual(q) for k, q in CON.items()}


def restrict(q, line):
    """Restrict a ternary form to the line a s + b t + c u = 0, as a binary form in (P, Q)."""
    A_, B_, C_ = line
    P, Q = sp.symbols('P Q')
    # a basis of the line: two independent solutions of A s + B t + C u = 0
    M = sp.Matrix([[A_, B_, C_]])
    ns = M.nullspace()
    v1, v2 = ns[0], ns[1]
    pt = [P * v1[i] + Q * v2[i] for i in range(3)]
    return sp.factor(sp.expand(q.subs({s: pt[0], t: pt[1], u: pt[2]}, simultaneous=True))), P, Q


if __name__ == '__main__':
    found = {}
    for pair in itertools.combinations('bdfh', 2):
        eqs = [DUAL[k] for k in pair]
        for chart in ({a: 1}, {a: 0, b_: 1}, {a: 0, b_: 0, c: 1}):
            e = [q.subs(chart) for q in eqs]
            free = [v for v in (a, b_, c) if v not in chart]
            try:
                sols = sp.solve(e, free, dict=True)
            except Exception:
                sols = []
            for so in sols:
                line = tuple(sp.nsimplify(chart.get(v, so.get(v, v))) for v in (a, b_, c))
                if any(x.free_symbols for x in line):
                    continue
                rational = all(x.is_rational for x in line)
                tang = [k for k in 'bdfh' if sp.simplify(DUAL[k].subs({a: line[0], b_: line[1], c: line[2]})) == 0]
                found[line] = (rational, tang)
    print(f'{len(found)} lines tangent to at least two conics (over Qbar)')
    for line, (rational, tang) in sorted(found.items(), key=lambda kv: (not kv[1][0], -len(kv[1][1]))):
        if not rational:
            print(f'  {str(line):<40} irrational; tangent to {tang}')
            continue
        info = []
        for k in 'bdfh':
            r, P, Q = restrict(CON[k], line)
            info.append(f'{k}|={r}')
        print(f'  {str(line):<40} RATIONAL; tangent to {tang}; ' + '; '.join(info))
