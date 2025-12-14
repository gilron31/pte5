import sympy as sp
import numpy as np
import math
import collections
import tqdm
import itertools
import random
from loguru import logger
import argparse
import cProfile
import pstats

import utils


class Enumerator:
    """
    - max_radius: The Radius is 2*P = a**2 + b**2.
    - min_primes: The Radius can be assumed to be composed of a minimal number of primes.
    """

    def __init__(self):
        self.skipped_norms = []
        self.parameterizers = {}

    def enumerate_radius_factorized(
        self,
        factorization,
        prime_decompositions,
        validate_2=True,
    ):
        if validate_2:
            assert 2 in factorization
            assert factorization[2] <= 2
        all_even_mults = True
        for p, m in factorization.items():
            if p == 2:
                continue
            assert p % 4 != 3
            all_even_mults &= m % 2 == 0
        if all_even_mults:
            self.skipped_norms.append(
                utils.calculate_norm_from_factorization(factorization)
            )
            return {}
        no_mults = all(v == 1 for v in factorization.values())
        num_primes = len(factorization)
        if no_mults:
            if num_primes not in self.parameterizers:
                self.parameterizers[num_primes] = (
                    utils.GaussianIntegersParameterization(num_primes)
                )
            points = list(
                set(
                    tuple(x)
                    for x in self.parameterizers[num_primes]
                    .get_pair_group_at(
                        [prime_decompositions[p] for p in factorization.keys()],
                        canonize_to_first_quadrant=True,
                    )
                    .transpose()
                    .tolist()
                )
            )
        else:
            points = utils.get_all_gaussian_integers_with_factored_norm(
                factorization,
                prime_decompositions={
                    p: prime_decompositions[p] for p in factorization.keys()
                },
            )
        Q_m_matches = utils.find_quadruplets(points)
        return Q_m_matches

    def enumerate_over_factorizations(
        self,
        prime_bound,
        num_primes,
        skip_squares=False,
        include_2=True,
        silently_square_everything=True,
        forced_primes=[],
    ):
        prime_decompositions = utils.get_all_decomposed_primes_up_to(prime_bound)
        if silently_square_everything:
            prime_decompositions = {
                p: (abs(d[0] ** 2 - d[1] ** 2), 2 * d[0] * d[1])
                for p, d in prime_decompositions.items()
            }
        primes = sorted(list(prime_decompositions.keys()))

        for p in forced_primes:
            primes.remove(p)

        radius_compositions = (
            itertools.combinations(primes[1:], num_primes)
            if skip_squares
            else itertools.combinations_with_replacement(primes[1:], num_primes)
        )

        results = {}
        for composition in tqdm.tqdm(radius_compositions):

            if include_2:
                composition = (2,) + composition
            composition = tuple(forced_primes) + composition
            factorization = dict(collections.Counter(composition))
            radius = utils.calculate_norm_from_factorization(factorization)
            matches = self.enumerate_radius_factorized(
                factorization,
                prime_decompositions,
                include_2,
            )
            if 134810 == radius:
                print(radius, factorization, matches)
            if len(matches) > 0:
                logger.info(f"{radius=} {factorization} {matches=}")
                for e4_sol in matches.values():
                    utils.analyze_E4_sol(e4_sol)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime_bound", type=int, required=True)
    parser.add_argument("--num_primes", type=int, required=True)
    parser.add_argument("--skip_squares", action="store_true")
    parser.add_argument("--include_2", action="store_true")
    parser.add_argument("--silently_square_everything", action="store_true")

    args = parser.parse_args()

    enum = Enumerator()

    # profiler = cProfile.Profile()

    # profiler.enable()
    enum.enumerate_over_factorizations(
        args.prime_bound,
        args.num_primes,
        skip_squares=args.skip_squares,
        include_2=args.include_2,
        silently_square_everything=args.silently_square_everything,
        forced_primes=[],
    )
    # profiler.disable()
    # stats = pstats.Stats(profiler)
    # stats.strip_dirs()
    # stats.sort_stats("cumulative")
    # stats.print_stats(100)
