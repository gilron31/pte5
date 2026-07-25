from code.enumeration_v1 import Enumerator
from loguru import logger


def test_perf_k4_varying_B():
    K = 4
    Bs = [200, 300, 400, 500]
    enum = Enumerator(max(Bs))
    for B in Bs:
        logger.info(f"{K=} {B=}")
        _, stats = enum.meet_points_from_factor_base_combinations(
            k=K, with_multiplicity=True, bound=B
        )
        logger.info(stats)


def test_perf_B200_varying_k():
    B = 200
    enum = Enumerator(B)
    for K in [3, 4, 5]:
        logger.info(f"{K=} {B=}")
        _, stats = enum.meet_points_from_factor_base_combinations(
            k=K,
            with_multiplicity=True,
        )
        logger.info(stats)
