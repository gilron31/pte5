from enumeration_v1 import Enumerator
import pytest


@pytest.mark.parametrize("with_multiplicity", [True, False])
def test_B200_k4(with_multiplicity):
    enum = Enumerator(200)
    res, _ = enum.meet_points_from_factor_base_combinations(
        k=4,
        with_multiplicity=with_multiplicity,
    )

    res = enum.enrich_results(res, print_analysis=False)
    assert len(res) == 2
    res = sorted(res, key=lambda x: x[0])
    assert sorted(res[0][0]) == [2, 5, 13, 17, 61]
    assert sorted(res[1][0]) == [2, 5, 29, 41, 101]
    # equating points is nabaz implement later


@pytest.mark.parametrize("with_multiplicity", [True, False])
def test_B200_k5(with_multiplicity):
    enum = Enumerator(200)
    res, _ = enum.meet_points_from_factor_base_combinations(
        k=5,
        with_multiplicity=with_multiplicity,
    )
    res = enum.enrich_results(res, print_analysis=False)
    assert len(res) == 3
    res = sorted(res, key=lambda x: x[0])
    assert sorted(res[0][0]) == [2, 5, 13, 17, 53, 101]
    assert sorted(res[1][0]) == [2, 5, 17, 41, 149, 157]
    assert sorted(res[2][0]) == [2, 13, 17, 61, 109, 113]
    # equating points is nabaz implement later
