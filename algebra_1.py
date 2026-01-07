import marimo

__generated_with = "0.18.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import sympy as sp
    return mo, sp


@app.cell
def _(sp):
    conj = lambda x: x.subs(sp.I, -sp.I)


    def as_complex(x):
        c = x.as_poly(sp.I).coeffs()
        return c[0] * sp.I + c[1]


    def as_tuple(x):
        c = x.as_poly(sp.I).coeffs()
        return c[0], c[1]


    def norm(z):
        return z[0] ** 2 + z[1] ** 2


    def calc_8L_2(p0, p1):
        A, B = as_tuple(p0)
        C, D = as_tuple(p1)
        return (A**2 - B**2) ** 2 + (C**2 - D**2) ** 2
    return as_tuple, calc_8L_2, conj, norm


@app.cell
def _(as_tuple, calc_8L_2, conj, mo, norm, sp):
    m, n, p, q = sp.var("m,n,p,q")
    z0 = m + sp.I * n
    z1 = p + sp.I * q
    AB = z0 * z1
    CD = z0 * conj(z1)
    Q = calc_8L_2(AB, CD)
    mo.vstack([norm(as_tuple(AB)), norm(as_tuple(CD)), Q.expand()])
    return AB, CD, Q


@app.cell
def _(as_tuple, calc_8L_2, conj, mo, norm, sp):
    m_, n_, p_, q_ = sp.var("m',n',p',q'")
    z0_ = m_ + sp.I * n_
    z1_ = p_ + sp.I * q_
    AB_ = z0_ * z1_
    CD_ = z0_ * conj(z1_)
    Q_ = calc_8L_2(AB_, CD_)
    mo.vstack(
        [norm(as_tuple(AB_)).expand(), norm(as_tuple(CD_)).expand(), Q_.expand()]
    )
    return AB_, CD_, Q_


@app.cell
def _(AB, AB_, CD, CD_, Q, Q_, as_tuple, mo, norm):
    Q_constr = Q_ - Q

    norm_constr = norm(as_tuple(AB)) - norm(as_tuple(AB_))

    constr_2_3 = (
        as_tuple(AB)[0]
        - as_tuple(CD)[0]
        - (as_tuple(AB_)[0] - as_tuple(CD_)[0])  # A - C = E - G
    )

    constr_2_4 = (  # A + B + C + D = E + F + G + H
        as_tuple(AB)[0]
        + as_tuple(AB)[1]
        + as_tuple(CD)[0]
        + as_tuple(CD)[1]
        - (
            as_tuple(AB_)[0]
            + as_tuple(AB_)[1]
            + as_tuple(CD_)[0]
            + as_tuple(CD_)[1]
        )
    )


    mo.vstack(
        [
            Q_constr.expand(),
            norm_constr.expand(),
            constr_2_3.simplify(),
            constr_2_4.simplify(),
        ]
    )
    return Q_constr, norm_constr


@app.cell
def _(Q_constr, norm_constr, sp):
    constraints = [
        Q_constr,
        norm_constr,
        # constr_2_3,
        # constr_2_4,
    ]

    r = sp.groebner(constraints)
    r.exprs
    return (r,)


@app.cell
def _(r, sp):
    f_r = [
        (x, sp.factor_list(x)) for x in r.exprs if len(sp.factor_list(x)[1]) > 1
    ]
    f_r
    return


@app.cell
def _():
    return


@app.cell
def _(mo):
    mo.md(r"""
    # Something inbetween me and Bremner
    """)
    return


@app.cell
def _(as_tuple, conj, mo, sp):
    alpha, beta, gamma, delta = sp.var(r"\alpha \beta,\gamma, \delta")
    A_, B_, C_, D_, E_, F_, G_, H_ = sp.var("A_,B_,C_,D_,E_,F_,G_,H_")
    P = sp.var("P")
    z0__ = alpha + sp.I * beta
    z1__ = gamma + sp.I * delta
    rhs = z0__ * z1__
    lhs = z0__ * conj(z1__)


    mo.vstack(
        [
            [
                A_ * B_ - as_tuple(rhs)[0],
                C_ * D_ - as_tuple(rhs)[1],
                E_ * F_ - as_tuple(lhs)[0],
                G_ * H_ - as_tuple(lhs)[1],
            ],
            [
                A_**2 - P - as_tuple(rhs)[0],
                C_**2 - P - as_tuple(rhs)[1],
                E_**2 - P - as_tuple(lhs)[0],
                G_**2 - P - as_tuple(lhs)[1],
            ],
            [
                A_**2 - B_**2 - as_tuple(rhs)[0],
                C_**2 - D_**2 - as_tuple(rhs)[1],
                E_**2 - F_**2 - as_tuple(lhs)[0],
                G_**2 - H_**2 - as_tuple(lhs)[1],
            ],
        ]
    )
    return A_, B_, C_, D_, E_, F_, G_, H_, alpha, beta, delta, gamma, lhs, rhs


@app.cell(hide_code=True)
def _(alpha, as_tuple, beta, delta, gamma, lhs, mo, rhs, sp):
    _s, _t, _u, _v = sp.var("s,t,u,v")
    subs_ = {alpha: _s + _t, beta: _u + _v, gamma: _s - _t, delta: _u - _v}
    mo.vstack([as_tuple(rhs.subs(subs_)), as_tuple(lhs.subs(subs_))])
    return (subs_,)


@app.cell
def _(alpha, as_tuple, beta, delta, gamma, lhs, mo, rhs):
    subs__ = {alpha: _s + _t, beta: _u + _v, gamma: _s - _t, delta: _u - _v}

    mo.vstack(
        [
            as_tuple(rhs)[0] - as_tuple(lhs)[0],
            as_tuple(rhs.subs(subs__))[0] - as_tuple(lhs.subs(subs__))[0],
        ]
    )
    return


@app.cell
def _(as_tuple, lhs, mo, rhs, subs_):
    mo.vstack(
        [
            as_tuple(rhs)[1] - as_tuple(lhs)[1],
            as_tuple(rhs.subs(subs_))[1] - as_tuple(lhs.subs(subs_))[1],
        ]
    )
    return


@app.cell
def _(C, D, G, H, as_tuple, lhs, mo, rhs, subs_):
    # sp.diophantine(AB -  as_tuple(lhs)[0])
    mo.vstack(
        [
            C**2 - D**2 - as_tuple(rhs.subs(subs_))[1],
            G**2 - H**2 - as_tuple(lhs.subs(subs_))[1],
        ]
    )
    return


@app.cell
def _(C, D, as_tuple, rhs, s, sp, subs_, t, u, v):
    _E1 = C**2 - D**2 - as_tuple(rhs.subs(subs_))[1]
    _x, _y, _z, _w, _e, _f, _k = sp.var("x, y, z, w, e, f, k")
    _subs = {
        C: 0 + _x * _k,
        t: 1 + _e * _k,
        u: 0 + _f * _k,
        D: 0 + _y * _k,
        s: 1 + _w * _k,
        v: 0 + _z * _k,
    }
    _E1_ = _E1.subs(_subs)
    _k_sol = sp.solve(_E1_, _k)[1]
    _k_denom = sp.denom(_k_sol)
    _final_param = {
        s: (v.subs({_k: _k_sol}) * _k_denom).simplify() for s, v in _subs.items()
    }
    assert _E1.subs(_final_param).simplify() == 0
    _final_param
    return


@app.cell
def _(C, D, G, H, as_tuple, lhs, rhs, s, sp, subs_, t):
    _E1 = C**2 - D**2 - as_tuple(rhs.subs(subs_))[1]
    _E2 = G**2 - H**2 - as_tuple(lhs.subs(subs_))[1]
    _E_sum = _E1 + _E2

    _x, _y, _z, _w, _e, _f, _k = sp.var("x, y, z, w, e, f, k")
    _subs = {
        C: 0 + _x * _k,
        t: 1 + _e * _k,
        G: 0 + _f * _k,
        D: 0 + _y * _k,
        s: 1 + _w * _k,
        H: 0 + _z * _k,
    }
    _E_sum_ = _E_sum.subs(_subs)
    _k_sol = sp.solve(_E_sum_, _k)[1]
    _k_denom = sp.denom(_k_sol)
    _final_param = {
        _s: (_v.subs({_k: _k_sol}) * _k_denom).simplify()
        for _s, _v in _subs.items()
    }
    assert _E_sum.subs(_final_param).simplify() == 0
    _final_param
    return


@app.cell
def _(A_, B_, C_, D_, E_, F_, G_, H_, as_tuple, lhs, rhs, sp):
    _res = sp.groebner(
        [
            A_ * B_ - as_tuple(rhs)[0],
            C_ * D_ - as_tuple(rhs)[1],
            E_ * F_ - as_tuple(lhs)[0],
            G_ * H_ - as_tuple(lhs)[1],
            (A_**2 + B_**2) - (C_**2 + D_**2),
            (A_**2 + B_**2) - (E_**2 + F_**2),
            (A_**2 + B_**2) - (G_**2 + H_**2),
        ]
    )
    # [sp.factor_list(x) for x in _res.exprs if len(sp.factor_list(x)[1]) > 1]
    _res.exprs
    return


@app.cell
def _(rhs):
    rhs.atoms()
    return


@app.cell
def _(
    A_,
    B_,
    C_,
    D_,
    E_,
    F_,
    G_,
    H_,
    alpha,
    as_tuple,
    beta,
    delta,
    gamma,
    lhs,
    mo,
    rhs,
    sp,
):
    _E1 = A_ * B_ - as_tuple(rhs)[0]

    _x, _y, _z, _w, _e, _f, _k = sp.var("x, y, z, w, e, f, k")
    _subs = {
        A_: 0 + _x * _k,
        B_: 0 + _e * _k,
        alpha: 1 + _f * _k,
        beta: 0 + _y * _k,
        gamma: -1 + _w * _k,
        delta: 0 + _z * _k,
    }
    _E1_ = _E1.subs(_subs)
    _k_sol = sp.solve(_E1_, _k)[1]
    _k_denom = sp.denom(_k_sol)
    _final_param = {
        s: (v.subs({_k: _k_sol}) * _k_denom).simplify() for s, v in _subs.items()
    }
    assert _E1.subs(_final_param).simplify() == 0
    # _final_param
    # _E1
    _some_res = [
        x.subs(_final_param).expand()
        for x in [
            A_ * B_ - as_tuple(rhs)[0],
            C_ * D_ - as_tuple(rhs)[1],
            E_ * F_ - as_tuple(lhs)[0],
            G_ * H_ - as_tuple(lhs)[1],
        ]
    ]

    mo.vstack([
        _final_param, _some_res
    ])
    return


@app.cell
def _(A_, B_, C_, D_, E_, F_, G_, H_, as_tuple, lhs, rhs):
    [
        x.subs(_final_param)
        for x in [
            A_ * B_ - as_tuple(rhs)[0],
            C_ * D_ - as_tuple(rhs)[1],
            E_ * F_ - as_tuple(lhs)[0],
            G_ * H_ - as_tuple(lhs)[1],
        ]
    ]
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Following up on Bremner's construction
    """)
    return


@app.cell
def _(sp):
    A, B, C, D, E, F, G, H = sp.var("A,B,C,D,E,F,G,H")
    s, t, u = sp.var("s,t,u")
    return A, B, C, D, E, F, G, H, s, t, u


@app.cell
def _(A, B, C, D, E, F, G, H, mo):
    constr_Q = (A * B) ** 2 + (C * D) ** 2 - ((E * F) ** 2 + (G * H) ** 2)
    constr_P_0 = (A**2 + B**2) - (C**2 + D**2)
    constr_P_1 = (A**2 + B**2) - (E**2 + F**2)
    constr_P_2 = (A**2 + B**2) - (G**2 + H**2)

    mo.vstack(
        [
            constr_P_0,
            constr_P_1,
            constr_P_2,
            constr_Q,
        ]
    )
    return constr_P_0, constr_P_1, constr_P_2, constr_Q


@app.cell
def _(A, C, E, G, constr_P_0, constr_P_1, constr_P_2, constr_Q, mo, s, t, u):
    eq_2_3_subs = {A: -s + t, C: s + t, E: -s + u, G: s + u}
    mo.vstack(
        [
            constr_P_0.subs(eq_2_3_subs).simplify(),
            constr_P_1.subs(eq_2_3_subs).expand(),
            constr_P_2.subs(eq_2_3_subs).expand(),
            constr_Q.subs(eq_2_3_subs).expand().simplify(),
        ]
    )
    return (eq_2_3_subs,)


@app.cell
def _(
    B,
    D,
    F,
    H,
    constr_P_0,
    constr_P_1,
    constr_P_2,
    constr_Q,
    eq_2_3_subs,
    sp,
):
    sols = sp.linsolve(
        (
            constr_P_0.subs(eq_2_3_subs).simplify(),
            constr_P_1.subs(eq_2_3_subs).expand(),
            constr_P_2.subs(eq_2_3_subs).expand(),
            constr_Q.subs(eq_2_3_subs).expand().simplify(),
        ),
        B**2,
        D**2,
        F**2,
        H**2,
    )
    sols
    return (sols,)


@app.cell
def _(constr_P_0, constr_P_1, constr_P_2, constr_Q, sp):
    _res = sp.groebner([constr_Q, constr_P_0, constr_P_1, constr_P_2])
    # sp.factor_list(_res.exprs[1])
    constr_Q_grob = _res.exprs[1]

    _res.exprs
    return (constr_Q_grob,)


@app.cell
def _(B, D, F, G, H, constr_Q_grob):
    _x = (
        (B**2 - G**2) * (B**2 - H**2)
        + (D**2 - G**2) * (D**2 - H**2)
        - (F**2 - G**2) * (F**2 - H**2)
    )
    assert (_x - constr_Q_grob).simplify() == 0
    _x.expand()
    return


@app.cell
def _(B, D, F, G, H):
    _y = B**4 + D**4 - F**4 + (G**2 + H**2) * (F**2 - B**2 - D**2) + G**2 * H**2
    _y.expand()
    return


@app.cell
def _(sp):
    x, y, z, w = sp.var("x,y,z,w")
    c_x, c_y, c_z, c_w = sp.var("c_x, c_y, c_z, c_w")
    return


@app.cell
def _(B, sols, sp):
    sp.diophantine(B**2 - list(sols)[0][0])
    # B**2 - list(sols)[0][0]
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
