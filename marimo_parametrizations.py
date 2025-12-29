import marimo

__generated_with = "0.18.4"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    import sympy as sp
    from collections import defaultdict
    import itertools

    return defaultdict, itertools, mo, sp


@app.cell
def _(mo, sp):
    N_POINTS = 8

    conj = lambda x: x.subs(sp.I, -sp.I)

    def as_complex(x):
        c = x.as_poly(sp.I).coeffs()
        return c[0] * sp.I + c[1]

    def as_tuple(x):
        c = x.as_poly(sp.I).coeffs()
        return c[0], c[1]

    ss = sp.var("".join([f"s{i}," for i in range(N_POINTS)]))
    ts = sp.var("".join([f"t{i}," for i in range(N_POINTS)]))

    checkboxes = mo.ui.array(
        [mo.ui.checkbox(label=f"{i}", value=False) for i in range(8)]
    )

    return as_tuple, checkboxes, conj, ss, ts


@app.cell
def _(checkboxes, mo):
    mo.vstack(checkboxes)
    return


@app.cell
def _(checkboxes, sp, ss, ts):
    ps = [
        s + sp.I * t if checkboxes[i].value else sp.Integer(1)
        for i, (s, t) in enumerate(zip(ss, ts))
    ]
    ps
    return (ps,)


@app.cell
def _(as_tuple, conj, ps):

    A0 = as_tuple(ps[0] * ps[1] * ps[2] * ps[3] * ps[4] * ps[5] * ps[6] * ps[7])
    A1 = as_tuple(
        ps[0]
        * conj(ps[1])
        * ps[2]
        * conj(ps[3])
        * ps[4]
        * conj(ps[5])
        * ps[6]
        * conj(ps[7])
    )
    A2 = as_tuple(
        ps[0]
        * ps[1]
        * conj(ps[2])
        * conj(ps[3])
        * ps[4]
        * ps[5]
        * conj(ps[6])
        * conj(ps[7])
    )
    A3 = as_tuple(
        ps[0]
        * ps[1]
        * ps[2]
        * ps[3]
        * conj(ps[4])
        * conj(ps[5])
        * conj(ps[6])
        * conj(ps[7])
    )

    Q12 = (A1[0] * A1[1]) ** 2 + (A2[0] * A2[1]) ** 2
    Q13 = (A1[0] * A1[1]) ** 2 + (A3[0] * A3[1]) ** 2
    Q10 = (A1[0] * A1[1]) ** 2 + (A0[0] * A0[1]) ** 2
    Q23 = (A2[0] * A2[1]) ** 2 + (A3[0] * A3[1]) ** 2
    Q20 = (A2[0] * A2[1]) ** 2 + (A0[0] * A0[1]) ** 2
    Q30 = (A3[0] * A3[1]) ** 2 + (A0[0] * A0[1]) ** 2

    D1 = (Q12 - Q30).as_poly()
    D2 = (Q13 - Q20).as_poly()
    D3 = (Q10 - Q23).as_poly()
    return A0, A1, A2, A3, D1


@app.cell
def _(checkboxes, mo):
    mo.vstack(checkboxes)
    return


@app.cell
def _(D1):
    D1

    return


@app.cell
def _(D1):
    D1.factor_list()
    return


@app.cell
def _(A0, A1, A2, A3, ss, ts):

    sub_dict = {ss[0]: 3, ss[7]: 3, ss[2]:4, ss[6]:6, ts[0]: -1, ts[7]:2, ts[2]:1, ts[6]:5}
    sol = [A0[0], A0[1], A1[0], A1[1], A2[0], A2[1], A3[0], A3[1]]
    [x.subs(sub_dict) for x in sol]
    return


@app.cell
def _(A0, ss, ts):
    A0[1].subs({ss[0]: 3, ss[7]: 3, ss[2]:4, ss[6]:6, ts[0]: -1, ts[7]:2, ts[2]:1, ts[6]:5})
    return


@app.cell
def _(A3, ss, ts):
    A3[0].subs({ss[0]: 3, ss[7]: 3, ss[2]:4, ss[6]:6, ts[0]: -1, ts[7]:2, ts[2]:1, ts[6]:5})
    return


@app.cell
def _(defaultdict, itertools):
    def get_all_values(pol, p, n_vars):
        rv = defaultdict(list)
        for v in itertools.product(range(p), repeat=n_vars):
            rv[pol(*v) % p].append(v)
        return rv

    def get_all_zeros(pol, p, n_vars):
        rv = []
        for v in itertools.product(range(p), repeat=n_vars):
            if pol(*v) % p == 0:
                rv.append(v)
        return rv
    return


@app.cell
def _():
    # rv = get_all_zeros(D1, 3, 8)
    # rv
    return


@app.cell
def _():
    3**8
    return


if __name__ == "__main__":
    app.run()
