from utils import verify_skew_normal_squareish
import pytest


def test_verify_skew_normal_squareish():
    assert verify_skew_normal_squareish([1, 2, 81])[1] == False
    assert verify_skew_normal_squareish([85, 4176, 2880**2])[1] == True
    assert (
        verify_skew_normal_squareish(
            [67405, 3525798096, 533470702551552000, 469208209191321600**2]
        )[1]
        == True
    )
