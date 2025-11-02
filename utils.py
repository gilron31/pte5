import sympy as sp
import numpy as np
import math
import collections
import tqdm
import itertools
import random
from loguru import logger


def compute_square_roots(l):
    rv = []
    for x in l:
        if x <= 0:
            return None
        r = sp.integer_nthroot(x, 2)
        if r[1]:
            rv.append(r[0])
        else:
            return None
    return rv


def verify_skew_normal_squareish(cs):
    iroots = [cs[-1]]
    iroots = compute_square_roots(iroots)

    is_gem = True
    for c in cs[-2::-1] + [0]:
        logger.debug(f"{c}, {iroots}")
        if iroots:
            iroots = [r for r_ in iroots for r in [c - r_, c + r_]]
        else:
            is_gem = False
            break
        if c != 0:
            iroots = compute_square_roots(iroots)
    logger.debug(f"{iroots=}")
    logger.info(f"{is_gem=}")
    return iroots, is_gem


def verify_skew_normal_symbolic(cs):
    x = sp.var("x")
    p = sp.poly(x)
    for c in cs:
        p = p**2 - c
    iroots = [r for r in p.all_roots() if r.is_Integer]
    is_gem = len(iroots) == 2 ** len(cs)
    print(f"{p=} \n{iroots=} \n{is_gem=}")
    return p, iroots, is_gem


class GaussianIntegersParameterization:
    def __init__(self, order):
        self.order = order
        self.guide = list(itertools.product([0, 1], repeat=self.order))
        self.re_guide = [g for g in self.guide if sum(g) % 2 == 0]
        self.im_guide = [g for g in self.guide if sum(g) % 2 == 1]
        self.im_sign_matrix = (
            -(
                (
                    (
                        (np.array(self.im_guide)[:, 1:])
                        @ np.array(
                            list(itertools.product([0, 1], repeat=order - 1))
                        ).transpose()
                    )
                    + ((np.array(self.im_guide).sum(axis=1, keepdims=True) - 1) // 2)
                )
                % 2
            )
            * 2
            + 1
        )
        self.re_sign_matrix = (
            -(
                (
                    (
                        (np.array(self.re_guide)[:, 1:])
                        @ np.array(
                            list(itertools.product([0, 1], repeat=order - 1))
                        ).transpose()
                    )
                    + ((np.array(self.re_guide).sum(axis=1, keepdims=True)) // 2)
                )
                % 2
            )
            * 2
            + 1
        )

    def get_pair_group_at(self, parameters, canonize_to_first_quadrant=False):
        assert len(parameters) == self.order
        re_monomials = [
            math.prod([p[g[i]] for i, p in enumerate(parameters)])
            for g in self.re_guide
        ]
        im_monomials = [
            math.prod([p[g[i]] for i, p in enumerate(parameters)])
            for g in self.im_guide
        ]

        im_coords = np.array(im_monomials) @ self.im_sign_matrix
        re_coords = np.array(re_monomials) @ self.re_sign_matrix

        rv = np.stack([re_coords, im_coords])

        if canonize_to_first_quadrant:
            return np.sort(np.abs(rv), axis=0)
        else:
            return rv


def decompose_prime(p):
    assert p % 4 == 1 or p == 2, "Prime is not decomposable"
    for x in range(int(math.sqrt(p)), 0, -1):
        y, valid = sp.integer_nthroot((p - x**2), 2)
        if valid:
            assert x >= y, (x, y)
            return (x, y)
    assert False, f"Could not decompose prime {p}"


def get_all_decomposable_primes_up_to(b):
    return [2] + [p for p in sp.primerange(b) if p % 4 == 1]


def get_all_decomposed_primes_up_to(b):
    primes = get_all_decomposable_primes_up_to(b)
    return {p: decompose_prime(p) for p in primes}


def get_all_gaussian_integers_with_norm(N):
    rv = []
    for a in range(math.floor(math.sqrt(N)) + 1):
        for b in range(
            math.floor(math.sqrt(N - a**2)), min(a, int(math.sqrt(N - a**2))) + 1
        ):
            if a**2 + b**2 == N:
                rv.append((a, b))
    return rv


def get_all_gaussian_integers_with_factored_norm(factors_list):
    gaussian_primes = []
    j_p_ranges = []
    one_mult_p_representative = None
    if 2 in factors_list:
        assert factors_list[2] == 1
    for p, e_p in factors_list.items():
        if p == 2:
            gaussian_primes.append((p, decompose_prime(p)))
            j_p_ranges.append([1])
            continue
        assert p % 4 == 1, f"{p}"
        if (one_mult_p_representative is None) and e_p % 2 == 1:
            one_mult_p_representative = p
            gaussian_primes.append((p, decompose_prime(p)))
            j_p_ranges.append(list(range(1, e_p + 1, 2)))
        else:
            gaussian_primes.append((p, decompose_prime(p)))
            j_p_ranges.append(list(range(-e_p, e_p + 1, 2)))

    assert one_mult_p_representative is not None
    rv = []

    for j_ps in itertools.product(*j_p_ranges):
        res_sp = 1
        for i, (p, d) in enumerate(gaussian_primes):
            d_sp = sp.ZZ_I(*d)
            d_conj_sp = sp.ZZ_I(*d[::-1])
            if j_ps[i] > 0:
                res_sp *= d_sp ** j_ps[i]
            else:
                res_sp *= d_conj_sp ** (-j_ps[i])
            res_sp *= p ** ((factors_list[p] - abs(j_ps[i])) // 2)
        rv.append(tuple(sorted([abs(res_sp.x), abs(res_sp.y)], reverse=True)))
    return sorted(rv)
