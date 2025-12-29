import marimo

__generated_with = "0.18.4"
app = marimo.App(width="full", layout_file="layouts/marimo_test.slides.json")


@app.cell
def _():
    import marimo as mo
    import utils
    import pandas as pd
    return mo, pd, utils


@app.cell
def _():
    return


@app.cell
def _(mo):
    B = mo.ui.number(step=1, value=101)
    B
    return (B,)


@app.cell
def _(B, pd, utils):
    fb = utils.get_all_decomposed_primes_up_to(B.value)
    df_fb = pd.DataFrame(fb)
    df_fb
    return (df_fb,)


@app.cell
def _(df_fb, mo):
    exp = mo.ui.dataframe(df_fb.transpose())
    exp


    return (exp,)


@app.cell
def _(exp):
    exp.value
    return


@app.cell
def _(mo):
    A = mo.ui.button(value=0, label=f"Press me", on_click=lambda x: x+1)
    A
    return (A,)


@app.cell
def _(A):
    A.value
    return


@app.cell
def _(pd):

    df = pd.DataFrame({"A": [1, 2, 3], "B": ["a", "b", "c"]})
    return (df,)


@app.cell
def _(df, mo):
    editor = mo.ui.data_editor(data=df, label="Edit Data")
    editor
    return


@app.cell
def _():
    return


@app.cell
def _(df):
    df
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
