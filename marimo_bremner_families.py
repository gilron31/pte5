import marimo

__generated_with = "0.18.4"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    import sympy as sp
    import utils
    from fractions import Fraction
    import pandas as pd
    import math
    return mo, pd, sp, utils


@app.cell
def _():
    return


@app.cell
def _():
    return


@app.cell
def _(mo):
    k_slider = mo.ui.number(value=5, step=1)
    k_slider
    return (k_slider,)


@app.cell
def _():
    return


@app.cell
def _(k_slider, utils):
    sol = utils.bremner_F2_point(k_slider.value)
    sol = utils.validate_and_canonize_E4_sol(sol)
    utils.analyze_E4_sol(sol, True)
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


@app.cell
def _(mo):

    k_slider_F1 = mo.ui.number(value=1, step=1)
    k_slider_F1
    return (k_slider_F1,)


@app.cell
def _(k_slider_F1, utils):

    _sol =  utils.bremner_F1_point( k_slider_F1.value)
    _sol = utils.validate_and_canonize_E4_sol(_sol)
    utils.analyze_E4_sol(_sol, True)
    return


@app.cell
def _(pd, utils):
    df = pd.DataFrame()
    # for _i in list(range(-13, -5)) + list(range(5, 14)):
    for k in range(5, 14):
        _sol = utils.bremner_F2_point(k)
        _sol = utils.validate_and_canonize_E4_sol(_sol)
        l = len(df)
        analysis = utils.analyze_E4_sol(_sol, factorize_gaussian_integers=True)
        df.loc[l, list(analysis.keys())] = list(analysis.values())
        df.loc[l, "A"] = _sol[0]
        df.loc[l, "B"] = _sol[1]
        df.loc[l, "C"] = _sol[2]
        df.loc[l, "D"] = _sol[3]
        df.loc[l, "E"] = _sol[4]
        df.loc[l, "F"] = _sol[5]
        df.loc[l, "G"] = _sol[6]
        df.loc[l, "H"] = _sol[7]
        df.loc[l, "Note"] = f"Bremner's Second family k={k}"

    for k in [1, 2, 3]:
        _sol = utils.bremner_F1_point(k)
        _sol = utils.validate_and_canonize_E4_sol(_sol)
        l = len(df)
        analysis = utils.analyze_E4_sol(_sol, factorize_gaussian_integers=True)
        df.loc[l, list(analysis.keys())] = list(analysis.values())
        df.loc[l, "A"] = _sol[0]
        df.loc[l, "B"] = _sol[1]
        df.loc[l, "C"] = _sol[2]
        df.loc[l, "D"] = _sol[3]
        df.loc[l, "E"] = _sol[4]
        df.loc[l, "F"] = _sol[5]
        df.loc[l, "G"] = _sol[6]
        df.loc[l, "H"] = _sol[7]
        df.loc[l, "Note"] = f"Bremner's First family k={k}"


    df.drop_duplicates(subset=["L_1"])
    return (df,)


@app.cell
def _(df):
    df[["L_1", "L_1_fact"]]
    return


@app.cell
def _(df):
    df[["Q_b", "Q_b_fact"]]
    return


@app.cell
def _(df, pd, sp):
    primes = set([k for x in df["Q_b_fact"] for k in x.keys() ])
    primes.remove(2)
    primes.remove(-1)
    primes

    primes_df = pd.DataFrame({"p" : list(primes)})
    primes_df["mod4"] = primes_df["p"] %4
    primes_df["mod8"] = primes_df["p"] %8
    primes_df["mod16"] = primes_df["p"] %16
    primes_df["factm1"] = primes_df["p"].apply(lambda p: sp.factorint(p-1))
    primes_df
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
