import marimo

__generated_with = "0.19.7"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    import code.utils as utils
    import sympy as sp
    import numpy as np
    import math
    import collections
    import tqdm
    import itertools
    import random
    from loguru import logger
    from fractions import Fraction

    return collections, random, sp, tqdm, utils


@app.cell
def _():
    return


@app.cell
def _(sp, utils):
    points = utils.comprime_first_octant_complex_up_to_norm(
        100000, exclude_even_norm=True
    )
    norms = [x.norm() for x in points]
    rpow4 = [(x**4).real for x in points]
    norm_factors = [sp.factorint(x) for x in norms]
    rpow4_factors = [sp.factorint(x) for x in rpow4]

    # def hash_point(point):
    #     norm = point.norm()
    #     rpow4 = (point**4).real
    #     norm_factors =
    len(points), len(points) ** 2 // 2
    return norms, points, rpow4


@app.cell
def _(random, sp):
    class AdditiveHash:
        def __init__(self, max_prime_hint=1000, seed=42, modulus=2**32):
            self.modulus = modulus
            self.rng = random.Random(seed)
            self.prime_weights = {}

            for i in sp.primerange(1, max_prime_hint + 1):
                self._assign_weight(i)

        def _assign_weight(self, p):
            if p not in self.prime_weights:
                curr_seed = self.rng.getrandbits(32) ^ p
                self.prime_weights[p] = curr_seed % self.modulus
            return self.prime_weights[p]

        def _get_weight(self, p):
            if p in self.prime_weights:
                return self.prime_weights[p]

            return hash(p) % self.modulus

        def hash(self, n):
            if n == 0:
                return 0
            if n == 1:
                return 0

            # sp.factorint returns a dictionary: {prime: multiplicity}
            # e.g., for 12, it returns {2: 2, 3: 1}
            factors = sp.factorint(n)

            total_hash = 0

            # Iterate over factors and sum their weights * exponents
            for p, exponent in factors.items():
                w = self._get_weight(p)
                term = (w * exponent) % self.modulus
                total_hash = (total_hash + term) % self.modulus

            return total_hash

    # --- Usage Example ---

    # Initialize with a hint to cover primes up to 100
    return (AdditiveHash,)


@app.cell
def _(AdditiveHash, collections, norms, points, rpow4, tqdm):
    def collision_to_quadruplet(i, j, points):
        g = points[i]
        h = points[j]
        return (g * h).to_first_octant(), (g * h.conjugate()).to_first_octant()

    def find_collisions(points, norms, rpow4):
        ht = collections.defaultdict(list)
        for i in tqdm.tqdm(range(len(points))):
            for j in range(i):
                norms_key = norms[i] * norms[j]
                rpow4_key = rpow4[i] * rpow4[j]
                key = (norms_key, rpow4_key)
                ht[key].append((i, j))
        return {k: v for k, v in ht.items() if len(v) > 1}

    def find_collisions_v2(points, norms, rpow4, hasher):
        ht0 = collections.defaultdict(list)
        hashs = [
            (hasher.hash(norms[i]) + hasher.hash(rpow4[i])) % hasher.modulus
            for i in range(len(points))
        ]
        print(f"first elimination")
        for i in tqdm.tqdm(range(len(points))):
            for j in range(i):
                key = (hashs[i] + hashs[j]) % hasher.modulus
                ht0[key].append((i, j))

        print(f"secondary elimination")
        ht0_candidates = {k: v for k, v in ht0.items() if len(v) > 1}
        ht1 = collections.defaultdict(list)
        for k, v in tqdm.tqdm(list(ht0_candidates.items())):
            for i, j in v:
                norms_key = norms[i] * norms[j]
                rpow4_key = rpow4[i] * rpow4[j]
                key = (norms_key, rpow4_key)
                ht1[key].append((i, j))
        return {k: v for k, v in ht1.items() if len(v) > 1}

    hasher = AdditiveHash(max_prime_hint=100, modulus=2**64)
    rs = find_collisions_v2(points, norms, rpow4, hasher)
    print(len(rs))
    rss = {
        k: (
            v,
            [(points[c[0]], points[c[1]]) for c in v],
            [collision_to_quadruplet(*c, points) for c in v],
        )
        for k, v in rs.items()
        if len(v) > 1
    }
    return rs, rss


@app.cell
def _():
    return


@app.cell
def _(rs):
    len(rs)
    return


@app.cell
def _(rs):
    {k: v for k, v in rs.items() if v == [(220, 88), (3945, 4)]}
    return


@app.cell
def _(rss):
    rss

    # next(iter(rs.items()))
    return


@app.cell
def _():
    [[220, 88], [3945, 4]]
    return


if __name__ == "__main__":
    app.run()
