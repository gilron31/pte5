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
    def __init__(
        self, bound, remove_2=True, debug=False, batch_size=4, hashtable_key_bits=15
    ):
        self.debug = debug
        self.bound = bound
        self.factor_base = {
            k: IntegerComplex(v[0], v[1])
            for k, v in get_all_decomposed_primes_up_to(self.bound, debug=debug).items()
            if (k != 2 or not remove_2)
        }
        self.max_factors = 10
        self.batch_size = batch_size
        self.points_scratch = np.zeros(
            (self.batch_size, 2, 2 ** (self.max_factors)), dtype=np.uint32
        )
        self.hashtable_key_size_bits = hashtable_key_bits
        self.hashtable: np.array = np.zeros(
            (self.batch_size, 2**self.hashtable_key_size_bits), dtype=np.uint8
        )
        self.indexing_dummy = np.arange(self.batch_size)[:, None]

    # @profile
    def get_all_points_from_factorization_approx(self, factors_batch):
        k = len(factors_batch[0])
        assert k <= self.max_factors

        initialized_factors = np.zeros((self.batch_size, 2, k))
        for i in range(len(factors_batch)):
            for j in range(k):
                initialized_factors[i, 0, j] = self.factor_base[
                    factors_batch[i][j]
                ].real
                initialized_factors[i, 1, j] = self.factor_base[
                    factors_batch[i][j]
                ].imag

        self.points_scratch[:, :, 1] = initialized_factors[:, :, 0]
        for i in range(k - 1):
            src_slice = self.points_scratch[:, :, 2 ** (i) : 2 ** (i + 1)]
            dst_slice_0 = self.points_scratch[
                :, :, 2 ** (i + 1) : 2 ** (i + 1) + 2 ** (i)
            ]
            dst_slice_1 = self.points_scratch[
                :, :, 2 ** (i + 1) + 2 ** (i) : 2 ** (i + 2)
            ]
            x = src_slice[:, 0, :] * initialized_factors[:, 0, i, None]
            y = src_slice[:, 0, :] * initialized_factors[:, 1, i, None]
            z = src_slice[:, 1, :] * initialized_factors[:, 0, i, None]
            w = src_slice[:, 1, :] * initialized_factors[:, 1, i, None]
            dst_slice_0[:, 0] = x - w
            dst_slice_0[:, 1] = y + z
            dst_slice_1[:, 0] = x + w
            dst_slice_1[:, 1] = y - z
        return self.points_scratch[:, :, 2 ** (k - 1) : 2 ** (k)]

    # @profile
    def meet_points_approx(self, points):
        point_projs = (points[0] * points[1]) ** 2

        # dedup
        self.hashtable[:] = 0
        for i in range(point_projs.shape[1]):
            for j in range(i):
                self.hashtable[
                    self.indexing_dummy,
                    (point_projs[:, i] + point_projs[:, j])
                    & ((2**self.hashtable_key_size_bits) - 1),
                ] += 1
        print(list(self.hashtable[0]))
        assert False
        return self.hashtable.max(axis=1) > 1
        # self.hashtable[self.indexing_dummy, point_projs] = 1

        # hashtable = dict()
        # for i, p in enumerate(point_projs):
        #     for j, q in enumerate(point_projs[:i]):
        #         v = p + q
        #         if v in hashtable:
        #             hashtable[v].append((i, j))
        #         else:
        #             hashtable[v] = [(i, j)]
        # return {k: v for k, v in hashtable.items() if len(v) > 1}

    # @profile
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

    # @profile
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

    # @profile
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

    def meet_points_from_factorization_hermetic(self, factors):
        points = self.get_all_points_from_factorization(factors)
        points_canonized = self.canonize_to_first_eighth(points)
        meet_result = self.meet_points(points_canonized)
        if len(meet_result) > 0:
            return factors, points_canonized, meet_result
        else:
            return None

    # @profile
    def meet_points_from_factorization(self, factors_batch):
        points_approx_batch = self.get_all_points_from_factorization_approx(
            factors_batch
        )

        meet_result_batch = self.meet_points_approx(points_approx_batch)

        print(sum(meet_result_batch))
        rv = []
        for i in range(len(factors_batch)):
            if meet_result_batch[i]:
                rv.append(
                    self.meet_points_from_factorization_hermetic(factors_batch[i])
                )
            else:
                rv.append(None)

        return rv

    @profile
    def meet_points_from_factorization_combinations(
        self, factors, k, with_multiplicity=True
    ):
        rv = []
        combinations = list(
            itertools.combinations_with_replacement(factors, k)
            if with_multiplicity
            else itertools.combinations(factors, k)
        )
        combination_count = len(combinations)

        start_time = time.time()

        num_batches = (combination_count + self.batch_size - 1) // self.batch_size

        with tqdm.tqdm(total=combination_count) as pbar:
            for batch_idx in range(num_batches):
                curr_batch_start = batch_idx * self.batch_size
                curr_batch_end = min(
                    combination_count, (batch_idx + 1) * self.batch_size
                )
                curr_batch_size = curr_batch_end - curr_batch_start
                choices = []
                for i in range(curr_batch_size):
                    choices.append(combinations[curr_batch_start + i])

                res = self.meet_points_from_factorization(choices)
                pbar.update(curr_batch_size)

                for res_ in res:
                    if res_:
                        rv.append(res_)

        elapsed_s = time.time() - start_time
        combinations_per_sec = combination_count / elapsed_s
        lg2_combinations_per_sec = math.log2(combinations_per_sec)
        # logger.debug(
        #     f"{lg2_combinations_per_sec=:0.2f}, {combination_count=} {elapsed_s=}"
        # )
        run_stats = {
            "combination_count": combination_count,
            "elapsed_s": elapsed_s,
            "combinations_per_sec": combinations_per_sec,
            "lg2_combinations_per_sec": lg2_combinations_per_sec,
        }

        return rv, run_stats

    def enrich_results(self, solutions, add_2=True, print_analysis=True):
        rv = []
        additional_factor = IntegerComplex(1, 1) if add_2 else IntegerComplex(1, 0)
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
                if print_analysis:
                    analyze_E4_sol(sol_points[-1])
            rv.append((factors, sol_points))
        return rv

    @profile
    def meet_points_from_factor_base_combinations(
        self, k, with_multiplicity=True, bound=None
    ):
        return self.meet_points_from_factorization_combinations(
            [p for p in self.factor_base.keys() if (p < bound if bound else True)],
            k,
            with_multiplicity=with_multiplicity,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime_bound", type=int, required=True)
    parser.add_argument("--num_primes", type=int, required=True)
    parser.add_argument("--with_multiplicity", action="store_true")

    args = parser.parse_args()

    enum = Enumerator(args.prime_bound, debug=True)

    res, stats = enum.meet_points_from_factor_base_combinations(
        args.num_primes,
        with_multiplicity=args.with_multiplicity,
    )
    logger.info(stats)
    res = enum.enrich_results(res)
    logger.info(res)
