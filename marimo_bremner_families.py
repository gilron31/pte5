import marimo

__generated_with = "0.18.4"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    import sympy as sp
    import utils
    from fractions import Fraction
    return Fraction, mo, sp, utils


@app.cell
def _():
    return


@app.cell
def _(Fraction, utils):
    coeffs = (Fraction(-1), Fraction(-8), Fraction(112))
    P = utils.GeneralPoint(Fraction(12), Fraction(40), coeffs)
    return (P,)


@app.cell
def _(mo):
    k_slider = mo.ui.number(value=1, step=1)
    k_slider
    return (k_slider,)


@app.cell
def _(P, k_slider):
    kP = P * k_slider.value
    kP
    return (kP,)


@app.cell
def _(kP, utils):
    sol = utils.point_to_bremner_2(kP)
    utils.analyze_E4_sol(
        (
            ((sol[0], sol[1]), (sol[2], sol[3])),
            ((sol[4], sol[5]), (sol[6], sol[7])),
        ),
        True,
    )
    return


@app.cell
def _(P, sp, utils):
    import pandas as pd

    df = pd.DataFrame(columns=["A", "B", "C", "D", "E", "F", "G", "H", "Note"])


    def validate_and_canonize_E4_sol(sol):
        sol = [abs(x) for x in sol]
        sol_pairs = [sorted(x, reverse=True) for x in zip(sol[::2], sol[1::2])]
        sol_quads = [
            sorted(x, reverse=True) for x in zip(sol_pairs[::2], sol_pairs[1::2])
        ]
        sol_octs = [
            sorted(x, reverse=True) for x in zip(sol_quads[::2], sol_quads[1::2])
        ]
        sol_quads = [x[i] for x in sol_octs for i in range(2)]
        sol_pairs = [x[i] for x in sol_quads for i in range(2)]
        sol = [x[i] for x in sol_pairs for i in range(2)]

        gcd = int(sp.gcd(sol))
        sol = [t // gcd for t in sol]
        radius = sol[0] ** 2 + sol[1] ** 2
        assert radius == sol[2] ** 2 + sol[3] ** 2
        assert radius == sol[4] ** 2 + sol[5] ** 2
        assert radius == sol[6] ** 2 + sol[7] ** 2
        Q_m = (sol[0] * sol[1]) ** 2 + (sol[2] * sol[3]) ** 2
        assert Q_m == (sol[4] * sol[5]) ** 2 + (sol[6] * sol[7]) ** 2
        if radius % 2 != 0:
            sol = [x * 2 for x in sol]
        return sol


    def analyze_E4_sol(
        sol, factorize_gaussian_integers=False
    ):
        radius = sol[0] ** 2 + sol[1] ** 2
        Q_m = (sol[0] * sol[1]) ** 2 + (sol[2] * sol[3]) ** 2
        Q_a = sol[0] ** 4 + sol[1] ** 4 + sol[2] ** 4 + sol[3] ** 4
        L_2 = (sol[0] ** 2 - sol[1] ** 2) ** 2 + (sol[2] ** 2 - sol[3] ** 2) ** 2
        Q_b = L_2 - 4 * Q_m

        assert L_2 % 8 == 0
        L_2 = L_2 // 8
        assert radius % 2 == 0
        L_1 = radius // 2

        rv = dict()

        rv["L_1"] = L_1
        rv["L_1_fact"] = sp.factorint(L_1)
        # rv["L_2"] = L_2
        # rv["L_2_fact"] = sp.factorint(L_2)
        # rv["Q_a"] = Q_a
        # rv["Q_a_fact"] = sp.factorint(Q_a)
        rv["Q_b"] = Q_b
        rv["Q_b_fact"] = sp.factorint(Q_b)
        # rv["Q_m"] = Q_m
        # rv["Q_m_fact"] = sp.factorint(Q_m)

        if factorize_gaussian_integers:
            i_part_00, factors_00, conj_mask_00 = utils.factorize_gaussian_integer(
                sol[0], sol[1]
            )
            i_part_01, factors_01, conj_mask_01 = utils.factorize_gaussian_integer(
                sol[2], sol[3]
            )
            i_part_10, factors_10, conj_mask_10 = utils.factorize_gaussian_integer(
                sol[4], sol[5]
            )
            i_part_11, factors_11, conj_mask_11 = utils.factorize_gaussian_integer(
                sol[6], sol[7]
            )

            assert i_part_00 == i_part_01
            assert i_part_00 == i_part_10
            assert i_part_00 == i_part_11

            assert factors_00 == factors_01
            assert factors_00 == factors_10
            assert factors_00 == factors_11

            xor_mask = lambda x, y: tuple(xx ^ yy for xx, yy in zip(x, y))
            xorred_00 = xor_mask(conj_mask_00, conj_mask_00)
            xorred_01 = xor_mask(conj_mask_00, conj_mask_01)
            xorred_10 = xor_mask(conj_mask_00, conj_mask_10)
            xorred_11 = xor_mask(conj_mask_00, conj_mask_11)
            syndrome_pattern = sorted(
                list(
                    set(
                        [
                            "".join([str(_x) for _x in x])
                            for x in zip(xorred_11, xorred_10, xorred_01)
                        ]
                    )
                )
            )

            rv["syndrome_pattern"] = ",".join(syndrome_pattern)
        return rv


    # for _i in list(range(-13, -5)) + list(range(5, 14)):
    for _i in range(5, 14):
        _sol = utils.point_to_bremner_2(P * _i)
        _sol = validate_and_canonize_E4_sol(_sol)
        analysis = analyze_E4_sol(_sol, factorize_gaussian_integers=True)
        df.loc[_i, list(analysis.keys())] = list(analysis.values())
        df.loc[_i, "A"] = _sol[0]
        df.loc[_i, "B"] = _sol[1]
        df.loc[_i, "C"] = _sol[2]
        df.loc[_i, "D"] = _sol[3]
        df.loc[_i, "E"] = _sol[4]
        df.loc[_i, "F"] = _sol[5]
        df.loc[_i, "G"] = _sol[6]
        df.loc[_i, "H"] = _sol[7]
        df.loc[_i, "Note"] = f"Bremner's Second family k={_i}"

    df
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
