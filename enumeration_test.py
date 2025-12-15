from enumeration_v1 import Enumerator
from utils import IntegerComplex


def test_basic_B200_k4():
    enum = Enumerator(200)
    res = enum.meet_points_from_factor_base(
        k=4,
        with_multiplicity=True,
    )

    res = enum.enrich_results(res)
    assert len(res) == 2
    res = sorted(res, key=lambda x: x[0])
    assert sorted(res[0][0]) == [2, 5, 13, 17, 61]
    assert sorted(res[1][0]) == [2, 5, 29, 41, 101]
    # equating points is nabaz implement later


def test_basic_B200_k5():
    enum = Enumerator(200)
    res = enum.meet_points_from_factor_base(
        k=5,
        with_multiplicity=True,
    )
    res = enum.enrich_results(res)
    assert len(res) == 3
    res = sorted(res, key=lambda x: x[0])
    assert sorted(res[0][0]) == [2, 5, 13, 17, 53, 101]
    assert sorted(res[1][0]) == [2, 5, 17, 41, 149, 157]
    assert sorted(res[2][0]) == [2, 13, 17, 61, 109, 113]
    # equating points is nabaz implement later
