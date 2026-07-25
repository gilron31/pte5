import marimo

__generated_with = "0.18.4"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    import sympy as sp
    import code.utils as utils
    from fractions import Fraction
    import pandas as pd
    import math
    import matplotlib.pyplot as plt
    import numpy as np
    from collections import defaultdict

    return defaultdict, mo, np, pd, plt, sp, utils


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

    df_bremner_sporadic = pd.DataFrame()
    sols = [
        [10043, 3167, 8473, 6253, 9677, 4153, 9217, 5093],
        [12233, 3661, 9689, 8317, 11549, 5447, 11287, 5971],
        [17923, 3295, 15277, 9935, 17573, 4825, 16025, 8677],
        [21523, 4567, 17593, 13213, 20687, 7493, 19483, 10223],
        [32507, 3031, 30109, 12623, 31799, 7397, 31037, 10129],
    ]
    for _i, _sol in enumerate(sols):
        _sol = utils.validate_and_canonize_E4_sol(_sol)
        l = len(df_bremner_sporadic)
        analysis = utils.analyze_E4_sol(
            _sol, factorize=False, factorize_gaussian_integers=False
        )
        df_bremner_sporadic.loc[l, list(analysis.keys())] = list(analysis.values())
        df_bremner_sporadic.loc[l, "Note"] = f"Bremner's sporadic i={_i}"

    df_gilro = pd.DataFrame()
    sols = [
        [61485, 249703, 148453, 209985, 167415, 195203, 46703, 252885],
        [3295, 17923, 9935, 15277, 4825, 17573, 8677, 16025],
        [7397, 31799, 10129, 31037, 3031, 32507, 12623, 30109],
        [8317, 9689, 3661, 12233, 5971, 11287, 5447, 11549],
        [1337, 3169, 1607, 3041, 2287, 2569, 809, 3343],
        [317, 1049, 541, 953, 719, 827, 139, 1087],
        [101, 353, 77, 359, 131, 343, 11, 367],
    ]
    for _i, _sol in enumerate(sols):
        _sol = utils.validate_and_canonize_E4_sol(_sol)
        l = len(df_gilro)
        analysis = utils.analyze_E4_sol(
            _sol, factorize=False, factorize_gaussian_integers=False
        )
        df_gilro.loc[l, list(analysis.keys())] = list(analysis.values())
        df_gilro.loc[l, "Note"] = f"gilro sporadic i={_i}"

    df = (
        pd.concat([df_bremner_first, df_bremner_second, df_bremner_sporadic, df_gilro])
        .reset_index(drop=True)
        .reset_index()
    )

    df = df.drop_duplicates(subset=["A", "B", "C", "D", "E", "F", "G", "H"])

    thresh = 2**200
    for _i, _row in df.iterrows():
        if _row["L_1"] < thresh:
            analysis = utils.analyze_E4_sol(
                _row[["A", "B", "C", "D", "E", "F", "G", "H"]].tolist(),
                factorize_gaussian_integers=_row["L_1"] < 2**120,
                factorize=True,
            )
            df.loc[_i, list(analysis.keys())] = list(analysis.values())

            # df.loc[_i, "syndrome_pattern"] = analysis["syndrome_pattern"]
            # df.at[_i, "L_1_fact"] = [analysis["L_1_fact"]]
            # df.at[_i, "Q_b_fact"] = [analysis["Q_b_fact"]]
    df[["Note", "L_1", "syndrome_pattern", "L_1_fact", "Q_b_fact"]]
    return df, df_bremner_first, df_bremner_second


@app.cell
def _(df):
    df["syndrome_pattern"].unique()
    # "".join([f"{k}^{v}\cdot " for k, v in df.loc[89]["Q_b_fact"].items()])
    return


app._unparsable_cell(
    r"""
    A	252885
    B	46703
    C	195203
    D	167415
    E	249703
    F	61485
    G	209985
    H	148453
    """,
    name="_",
)


@app.cell
def _(defaultdict):
    R = 1000
    invert = defaultdict(list)

    for i in range(1, R):
        for j in range(1 + i % 2, i, 2):
            # if sp.gcd(i, j)!= 1:
            #     continue
            rad = i**2 + j**2
            Q_b = i**4 + j**4 - 6 * (i**2) * (j**2)
            invert[(Q_b)].append((i, j))
    return (invert,)


@app.cell
def _(invert, sp):
    [
        (k, sp.factorint(k), v)
        # (k, v)
        for k, v in list(invert.items())[-100:]
        # if len(v) > 1 and all(sp.gcd(i, j) == 1 for i, j in v)
    ]
    return


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
def _(mo):
    mod_chooser = mo.ui.number(value=2, step=1)
    mod_chooser
    return (mod_chooser,)


@app.cell
def _(df, mod_chooser, pd, sp):
    Q_b_fact_primes = set([k for x in df["Q_b_fact"].dropna() for k in x.keys()])
    Q_b_fact_primes.remove(2)
    Q_b_fact_primes.remove(-1)
    Q_b_fact_primes

    primes_df = pd.DataFrame({"p": list(Q_b_fact_primes)})
    primes_df["mod4"] = primes_df["p"] % 4
    primes_df["mod8"] = primes_df["p"] % 8
    primes_df["mod"] = primes_df["p"] % mod_chooser.value
    primes_df["factm1"] = primes_df["p"].apply(lambda p: sp.factorint(p - 1))
    primes_df
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
