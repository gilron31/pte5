import sympy as sp
import numpy as np
import math
import collections
import tqdm
import itertools
import random
from loguru import logger

import utils


class Enumerator:
    """
    - max_radius: The Radius is 2*P = a**2 + b**2.
    - min_primes: The Radius can be assumed to be composed of a minimal number of primes.
    """

    def __init__(self):
        self.skipped_norms = []

    def enumerate_radius_factorized(self, factorization):
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
        points = utils.get_all_gaussian_integers_with_factored_norm(factorization)
        Q_m_matches = utils.find_quadruplets(points)
        return Q_m_matches

    def enumerate_over_factorizations(self, prime_bound, num_primes, skip_squares = False):
        prime_decompositions = utils.get_all_decomposed_primes_up_to(prime_bound)
        primes = sorted(list(prime_decompositions.keys()))
        
        if skip_squares:
        radius_factorizations = [
            dict(collections.Counter(x + (2,)))
            for x in itertools.combinations_with_replacement(primes[1:], num_primes)
        ]
        results = {}
        for factorization in tqdm.tqdm(radius_factorizations):
            radius = utils.calculate_norm_from_factorization(factorization)
            matches = self.enumerate_radius_factorized(factorization)
            if 134810 == radius:
                print(radius, factorization, matches)
            if len(matches) > 0:
                logger.info(f"{radius=} {factorization} {matches=}")
            # print(self.skipped_norms)


if __name__ == "__main__":
    enum = Enumerator()
    # enum.enumerate_radius_factorized(sp.factorint(2 * 5 * 13 * 17 * 29 * 37 * 53))
    enum.enumerate_over_factorizations(101, 8)
