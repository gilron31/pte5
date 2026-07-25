import marimo

__generated_with = "0.18.4"
app = marimo.App()


@app.cell
def _():
    import code.utils as utils
    import sympy as sp
    import itertools
    from collections import defaultdict
    import marimo as mo
    import numpy as np
    import pandas as pd
    import altair as alt

    return alt, defaultdict, itertools, mo, pd, sp


@app.cell
def _():
    # fb = utils.get_all_decomposed_primes_up_to(1000)
    # list(fb.keys())
    return


@app.cell
def _(mo, sp):
    x, y = sp.var("x, y")
    pol = sp.poly(x**2)
    p_slider = mo.ui.number(step=1, value=7)
    p_slider
    return p_slider, x, y


@app.cell
def _(defaultdict, itertools, p_slider, sp, x, y):
    p = p_slider.value

    def get_all_values(pol, p, n_vars):
        rv = defaultdict(list)
        for v in itertools.product(range(p), repeat=n_vars):
            rv[pol(*v) % p].append(v)
        return rv

    def get_all_values_unordered(pol, p, n_vars):
        rv = defaultdict(set)
        for v in itertools.product(range(p), repeat=n_vars):
            rv[pol(*v) % p].add(tuple(sorted(v, reverse=True)))
        return rv

    rv = get_all_values_unordered(sp.poly(x**2 + y**2), p, 2)
    return get_all_values_unordered, rv


@app.cell
def _(pd, rv):
    df = pd.DataFrame(rv[0], columns=["x", "y"]).sort_values(by="x", ignore_index=True)
    df
    return (df,)


@app.cell
def _(mo):
    mo.md(r"""
    ### Decomposable primes exhibit Criss-Cross Pattern

    Different primes show different slopes...
    """)
    return


@app.cell
def _(alt, df, mo):
    chart = mo.ui.altair_chart(alt.Chart(df).mark_point().encode(x="x", y="y"))
    chart
    return (chart,)


@app.cell
def _(chart):
    chart.value
    return


@app.cell
def _(defaultdict, get_all_values_unordered, sp, x, y):

    def Q_experiment_1(p_):
        points_by_norm_value = get_all_values_unordered(sp.poly(x**2 + y**2), p_, 2)

        # Iterate only over point pairs with the same norm
        def point_evaluation(x, y):
            return (x * y) ** 2 % p_

        rv = defaultdict(set)
        for norm_value, points in points_by_norm_value.items():
            points_l = list(points)
            for i, point_a in enumerate(points_l):
                a_eval = point_evaluation(*point_a)
                for j, point_b in enumerate(points_l[:i]):
                    b_eval = point_evaluation(*point_b)
                    if (a_eval + b_eval) % p_ == 0:
                        rv[norm_value].add(
                            tuple(sorted((point_a, point_b), key=lambda x: x[0]))
                        )
        return rv

    # This is a strict requirement since we require both P and Q_m to be divisible by the same prime
    def get_all_Q_zeros(p):
        point_with_zero_norm = get_all_values_unordered(sp.poly(x**2 + y**2), p, 2)[0]

        def point_evaluation(x, y):
            return (x * y) ** 2 % p

        eval_sum_values = defaultdict(set)
        points_l = list(point_with_zero_norm)
        for i, point_a in enumerate(points_l):
            a_eval = point_evaluation(*point_a)
            for j, point_b in enumerate(points_l[:i]):
                b_eval = point_evaluation(*point_b)
                eval_sum_values[(a_eval + b_eval) % p].add(
                    tuple(sorted((point_a, point_b), key=lambda x: x[0]))
                )
        return eval_sum_values

    return (Q_experiment_1,)


@app.cell
def _(mo, sp):
    primes = list(sp.primerange(1000))
    prime_selector = mo.ui.number(step=1, start=1)
    prime_selector
    return prime_selector, primes


@app.cell
def _(Q_experiment_1, prime_selector, primes):
    prime = primes[prime_selector.value]
    # prime = 233
    print(prime, prime % 4 == 1)
    Q_experiment_1(prime)
    return


@app.cell
def _(mo):
    mo.md(r"""
    # If we want P and Q_m to be divisible by a common factor

    Than than factor seems to have to be 1 mod 8 (not 4 which is trivial)
    """)
    return


@app.cell
def _(Q_experiment_1, primes):
    l = []
    num_primes = 63
    for p_ in primes[:num_primes]:
        res = Q_experiment_1(p_)
        if p_ % 8 == 1:
            assert len(res[0]) > 0, p_
        if len(res[0]) > 0:
            l.append(p_)

    l, [p_ for p_ in primes[:num_primes] if p_ % 8 == 1]
    return


@app.cell
def _(primes):
    primes[62]
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
