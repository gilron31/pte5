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
from line_profiler import profile

from utils import IntegerComplex, get_all_decomposed_primes_up_to, analyze_E4_sol


class Enumerator:
    def __init__(self, bound, remove_2=True):
        self.bound = bound
        self.factor_base = {
            k: IntegerComplex(v[0], v[1])
            for k, v in get_all_decomposed_primes_up_to(self.bound).items()
            if (k != 2 or not remove_2)
        }

    @profile
    def get_all_points_from_factorization(self, factors):
        n_factors = len(factors)
        all_points = []
        for i in range(n_factors):
            g = self.factor_base[factors[i]]
            # g_conj = g.conjugate()
            if len(all_points) > 0:
                # A trick for saving multiplications
                all_points_ = []
                for p in all_points:
                    x = p.real * g.real
                    y = p.real * g.imag
                    z = p.imag * g.real
                    w = p.imag * g.imag
                    all_points_.append(IntegerComplex(x - w, y + z))
                    all_points_.append(IntegerComplex(x + w, y - z))
                all_points = all_points_

                # The naive way:

                # all_points = [g * p for p in all_points] + [
                #     g_conj * p for p in all_points
                # ]
            else:
                all_points.append(g)
        return all_points

    @profile
    def canonize_to_first_eighth(self, points, dedup=True):
        rv = set() if dedup else list()
        for p in points:
            x = abs(p.real)
            y = abs(p.imag)
            if dedup:
                rv.add((min(x, y), max(x, y)))
            else:
                rv.append((min(x, y), max(x, y)))
        return [IntegerComplex(v[0], v[1]) for v in rv]

    @profile
    def meet_points(self, points):
        point_projs = [(p.real * p.imag) ** 2 for p in points]
        hashtable = dict()
        for i, p in enumerate(point_projs):
            for j, q in enumerate(point_projs[:i]):
                v = p + q
                if v in hashtable:
                    hashtable[v].append((i, j))
                else:
                    hashtable[v] = [(i, j)]
        return {k: v for k, v in hashtable.items() if len(v) > 1}

    @profile
    def meet_points_from_factorization(self, factors, k=None, with_multiplicity=True):
        if k is None:
            points = self.get_all_points_from_factorization(factors)
            points = self.canonize_to_first_eighth(points)
            return (
                factors,
                points,
                self.meet_points(points),
            )
        else:
            rv = []
            combinations = (
                itertools.combinations_with_replacement(factors, k)
                if with_multiplicity
                else itertools.combinations(factors, k)
            )
            for choice in tqdm.tqdm(list(combinations)):
                res = self.meet_points_from_factorization(choice)
                if len(res[2]) > 0:
                    rv.append(res)
            return rv

    @profile
    def enrich_results(self, solutions, add_2=True):
        rv = []
        additional_factor = IntegerComplex(1, 1) if add_2 else IntegerComplex(1, 0)
        # additional_factor = IntegerComplex(1, 0)
        for factors, points, meets in solutions:
            if add_2:
                factors = [2] + list(factors)
            sol_points = []
            for _, meet in meets.items():
                sol_points_ = self.canonize_to_first_eighth(
                    [
                        additional_factor * points[meet[0][0]],
                        additional_factor * points[meet[0][1]],
                        additional_factor * points[meet[1][0]],
                        additional_factor * points[meet[1][1]],
                    ],
                    dedup=False,
                )
                sol_points.append(
                    (
                        (sol_points_[0], sol_points_[1]),
                        (sol_points_[2], sol_points_[3]),
                    )
                )
                analyze_E4_sol(sol_points[-1])
            rv.append((factors, sol_points))
        return rv

    @profile
    def meet_points_from_factor_base(self, k, with_multiplicity=True):
        return self.meet_points_from_factorization(
            list(self.factor_base.keys()), k, with_multiplicity=with_multiplicity
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime_bound", type=int, required=True)
    parser.add_argument("--num_primes", type=int, required=True)
    parser.add_argument("--with_multiplicity", action="store_true")

    args = parser.parse_args()

    enum = Enumerator(args.prime_bound)

    res = enum.meet_points_from_factor_base(
        args.num_primes,
        with_multiplicity=args.with_multiplicity,
    )
    res = enum.enrich_results(res)
    logger.info(res)
