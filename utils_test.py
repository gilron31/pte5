from utils import (
    verify_skew_normal_squareish,
    verify_skew_normal_symbolic,
    get_all_gaussian_integers_with_norm,
    get_all_decomposed_primes_up_to,
    get_all_gaussian_integers_with_factored_norm,
)
import pytest
import sympy as sp


def test_verify_skew_normal_squareish():
    iroots, is_gem = verify_skew_normal_squareish([1, 2, 3])
    assert is_gem == False
    assert iroots == None
    iroots, is_gem = verify_skew_normal_squareish([85, 4176, 2880**2])
    assert is_gem == True
    assert sorted(iroots) == [-13, -11, -7, -1, 1, 7, 11, 13]
    iroots, is_gem = verify_skew_normal_squareish(
        [67405, 3525798096, 533470702551552000, 469208209191321600**2]
    )
    assert sorted(iroots) == [
        -367,
        -359,
        -353,
        -343,
        -131,
        -101,
        -77,
        -11,
        11,
        77,
        101,
        131,
        343,
        353,
        359,
        367,
    ]
    assert is_gem == True


def test_verify_skew_normal_symbolic():
    _, iroots, is_gem = verify_skew_normal_symbolic([1, 2, 3])
    assert is_gem == False
    assert iroots == []
    _, iroots, is_gem = verify_skew_normal_symbolic([85, 4176, 2880**2])
    assert is_gem == True
    assert sorted(iroots) == [-13, -11, -7, -1, 1, 7, 11, 13]
    _, iroots, is_gem = verify_skew_normal_symbolic(
        [67405, 3525798096, 533470702551552000, 469208209191321600**2]
    )
    assert sorted(iroots) == [
        -367,
        -359,
        -353,
        -343,
        -131,
        -101,
        -77,
        -11,
        11,
        77,
        101,
        131,
        343,
        353,
        359,
        367,
    ]
    assert is_gem == True


def test_decomposed_primes():
    B = 1000
    ps_and_decompositions = get_all_decomposed_primes_up_to(B)
    for p, d in ps_and_decompositions.items():
        exp_d = get_all_gaussian_integers_with_norm(p)
        assert len(exp_d) == 1, p
        assert exp_d[0] == d, p


def test_get_all_gaussian_integers_with_norm():
    N = 5
    expected = get_all_gaussian_integers_with_norm(N)
    result = get_all_gaussian_integers_with_factored_norm(sp.factorint(N))
    assert expected == result

    N = 2 * 5 * 13 * 17 * 29
    expected = get_all_gaussian_integers_with_norm(N)
    result = get_all_gaussian_integers_with_factored_norm(sp.factorint(N))
    assert expected == result

    N = 2 * 5 * 17 * 17 * 17 * 29 * 29
    expected = get_all_gaussian_integers_with_norm(N)
    result = get_all_gaussian_integers_with_factored_norm(sp.factorint(N))
    assert expected == result
