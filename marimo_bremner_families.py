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
    import matplotlib.pyplot as plt
    import numpy as np
    return mo, np, pd, plt, sp, utils


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
    _sol = utils.bremner_F1_point(k_slider_F1.value)
    _sol = utils.validate_and_canonize_E4_sol(_sol)
    utils.analyze_E4_sol(_sol, True)
    return


@app.cell
def _(pd, utils):
    df_bremner_first = pd.DataFrame()
    for k in range(1, 30):
        _sol = utils.bremner_F1_point(k)

        _sol = utils.validate_and_canonize_E4_sol(_sol)
        l = len(df_bremner_first)
        analysis = utils.analyze_E4_sol(
            _sol, factorize=False, factorize_gaussian_integers=False
        )
        df_bremner_first.loc[l, list(analysis.keys())] = list(analysis.values())
        df_bremner_first.loc[l, "Note"] = f"Bremner's Second family k={k}"

    df_bremner_second = pd.DataFrame()
    for k in range(5, 60):
        _sol = utils.bremner_F2_point(k)
        _sol = utils.validate_and_canonize_E4_sol(_sol)
        l = len(df_bremner_second)
        analysis = utils.analyze_E4_sol(
            _sol, factorize=False, factorize_gaussian_integers=False
        )
        df_bremner_second.loc[l, list(analysis.keys())] = list(analysis.values())
        df_bremner_second.loc[l, "Note"] = f"Bremner's First family k={k}"


    df = pd.concat([df_bremner_first, df_bremner_second]).reset_index(drop=True).reset_index()

    df = df.drop_duplicates(subset=["A", "B", "C", "D", "E", "F", "G", "H"])
    df
    return df, df_bremner_first, df_bremner_second


@app.cell
def _(df):
    df["L_1"].is_unique, df["Q_b"].is_unique
    return


@app.cell
def _(df_bremner_first):
    df_bremner_first
    return


@app.cell
def _(df_bremner_first, np, plt):
    _x = np.log2(df_bremner_first.index + 1e-3)[20:]
    _y = np.log2(df_bremner_first["lg2_coeff_size"].astype(int))[20:]
    print(np.polyfit(_x, _y, 1, full=True))
    plt.plot(_x, _y, "-*")


    return


@app.cell
def _(df_bremner_second, np, plt):

    _x = np.log2(df_bremner_second.index + 1e-3)[30:]
    _y = np.log2(df_bremner_second["lg2_coeff_size"].astype(int))[30:]
    print(np.polyfit(_x, _y, 1, full=True))
    plt.plot(_x, _y, "-*")
    return


@app.cell
def _(df):
    df["L_1"].is_unique, df["Q_b"].is_unique
    return


@app.cell
def _():
    return


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
    Q_b_fact_primes = set([k for x in df["Q_b_fact"] for k in x.keys()])
    Q_b_fact_primes.remove(2)
    Q_b_fact_primes.remove(-1)
    Q_b_fact_primes

    primes_df = pd.DataFrame({"p": list(Q_b_fact_primes)})
    primes_df["mod4"] = primes_df["p"] % 4
    primes_df["mod8"] = primes_df["p"] % 8
    primes_df["mod16"] = primes_df["p"] % 16
    primes_df["factm1"] = primes_df["p"].apply(lambda p: sp.factorint(p - 1))
    primes_df
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
