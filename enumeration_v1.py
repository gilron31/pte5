import sympy as sp
import numpy as np
import math
import collections
import tqdm
import itertools
import random
from loguru import logger
import argparse
import time
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
        self.max_factors = 10
        self.points_scratch = np.zeros((2, 2 ** (self.max_factors)), dtype=np.uint32)

    @profile
    def get_all_points_from_factorization_approx(self, factors):
        n_factors = len(factors)
        assert n_factors <= self.max_factors
        self.points_scratch[0, 1] = self.factor_base[factors[0]].real
        self.points_scratch[1, 1] = self.factor_base[factors[0]].imag
        for i in range(n_factors - 1):
            src_slice = self.points_scratch[:, 2 ** (i) : 2 ** (i + 1)]
            dst_slice_0 = self.points_scratch[:, 2 ** (i + 1) : 2 ** (i + 1) + 2 ** (i)]
            dst_slice_1 = self.points_scratch[:, 2 ** (i + 1) + 2 ** (i) : 2 ** (i + 2)]
            x = src_slice[0] * self.factor_base[factors[i + 1]].real
            y = src_slice[0] * self.factor_base[factors[i + 1]].imag
            z = src_slice[1] * self.factor_base[factors[i + 1]].real
            w = src_slice[1] * self.factor_base[factors[i + 1]].imag
            dst_slice_0[0] = x - w
            dst_slice_0[1] = y + z
            dst_slice_1[0] = x + w
            dst_slice_1[1] = y - z
        return self.points_scratch[:, 2 ** (n_factors - 1) :]

    @profile
    def meet_points_approx(self, points, dedup=True):
        point_projs = (points[0] * points[1]) ** 2

        if dedup:
            point_projs = list(set(point_projs))
        else:
            point_projs = list(point_projs)

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
            points_approx = self.get_all_points_from_factorization_approx(factors)
            meet_result = self.meet_points_approx(points_approx)
            if len(meet_result) == 0:
                return None

            points = self.get_all_points_from_factorization(factors)
            points = self.canonize_to_first_eighth(points)
            meet_result = self.meet_points(points)
            if len(meet_result) == 0:
                return None
            else:

                return (
                    factors,
                    points,
                    meet_result,
                )
        else:
            rv = []
            combinations = list(
                itertools.combinations_with_replacement(factors, k)
                if with_multiplicity
                else itertools.combinations(factors, k)
            )
            combination_count = len(combinations)
            start_time = time.time()
            for choice in tqdm.tqdm(combinations):
                res = self.meet_points_from_factorization(choice)
                if res:
                    rv.append(res)
            elapsed_s = time.time() - start_time
            combinations_per_sec = combination_count / elapsed_s
            lg2_combinations_per_sec = math.log2(combinations_per_sec)
            logger.info(
                f"{lg2_combinations_per_sec=:0.2f}, {combination_count=} {elapsed_s=}"
            )
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
