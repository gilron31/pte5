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
