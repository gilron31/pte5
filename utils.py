import sympy as sp
import numpy as np
import math
import collections
import tqdm
import itertools
import random
from loguru import logger
from fractions import Fraction


class IntegerComplex:
    def __init__(self, real, imag):
        self.real = real
        self.imag = imag

    def __repr__(self):
        return f"({self.real}, {self.imag})"

    def __str__(self):
        if self.imag >= 0:
            return f"{self.real} + {self.imag}i"
        else:
            # Handle negative imaginary part cleanly
            return f"{self.real} - {-self.imag}i"

    def __add__(self, other):
        return IntegerComplex(self.real + other.real, self.imag + other.imag)

    def __sub__(self, other):
        return IntegerComplex(self.real - other.real, self.imag - other.imag)

    def __mul__(self, other):
        if isinstance(other, IntegerComplex):
            new_real = (self.real * other.real) - (self.imag * other.imag)
            new_imag = (self.real * other.imag) + (self.imag * other.real)
            return IntegerComplex(new_real, new_imag)
        elif isinstance(other, int):
            return IntegerComplex(other * self.real, other * self.imag)

    def __pow__(self, pow):
        assert isinstance(pow, int)
        if pow == 0:
            return IntegerComplex(1, 0)
        assert pow > 0

        rv = self
        for _ in range(pow - 1):
            rv *= self
        return rv

    def __getitem__(self, index):
        if index == 0:
            return self.real
        if index == 1:
            return self.imag
        assert False, "Index out of range."

    def to_tuple(self):
        return (self.real, self.imag)

    def __hash__(self):
        return hash(self.to_tuple())

    def __eq__(self, other):
        assert isinstance(other, IntegerComplex)
        return self.real == other.real and self.imag == other.imag

    def __neg__(self):
        return IntegerComplex(-self.real, -self.imag)

    def conjugate(self):
        return IntegerComplex(self.real, -self.imag)

    def norm(self):
        return self.real**2 + self.imag**2


class GeneralPoint:
    def __init__(self, x, y, a_coeffs):
        """a_coeffs = (a2, a4, a6) for curve: y^2 = x^3 + a2*x^2 + a4*x + a6"""
        self.x = x
        self.y = y
        self.a2, self.a4, self.a6 = a_coeffs
        assert (
            y**2 - x**3 - self.a2 * x**2 - self.a4 * x - self.a6 == 0
        ), f"{y**2=},  {x**3 + self.a2 * x**2 + self.a4 * x + self.a6=}"

    def __add__(self, other):
        if self.x is None:
            return other
        if other.x is None:
            return self

        # Point negation for y^2 = f(x) is (x, -y)
        if self.x == other.x and self.y == -other.y:
            return GeneralPoint(None, None, (self.a2, self.a4, self.a6))

        if self.x != other.x:
            l = (other.y - self.y) / (other.x - self.x)
        else:
            # Doubling: Derivative includes the 2*a2*x term
            l = (3 * self.x**2 + 2 * self.a2 * self.x + self.a4) / (2 * self.y)

        # Main difference: subtract a2 when calculating x3
        x3 = l**2 - self.a2 - self.x - other.x
        y3 = l * (self.x - x3) - self.y

        return GeneralPoint(x3, y3, (self.a2, self.a4, self.a6))

    def __mul__(self, k):
        assert isinstance(k, int)
        assert k != 0, "Not supported"
        if k < 0:
            return GeneralPoint(self.x, -self.y, (self.a2, self.a4, self.a6)) * -k
        rv = self
        for i in range(k - 1):
            rv += self
        return rv

    def __repr__(self):
        return f"({self.x}, {self.y})"


"""
Magma in:

`http://magma.maths.usyd.edu.au/calc/`

```
PP<x,y,z>:=ProjectiveSpace(Rationals(),2);
C:=Curve(PP,3*x^2*y-4*x*y^2+y^3+(x^2 -14*x*y +9*y^2)*z +10*(x+y)*z^2);
P0:=C![1,1,0];
E, phi:=EllipticCurve(C,P0);
Em, psi:= MinimalModel(E);
E;
Em;
phi;
psi;
Rank(Em);

Inverse(psi);
is_bir, phi_inv := IsInvertible(phi);
phi_inv;
```

Magma out

```
Elliptic Curve defined by y^2 = x^3 + 4/5*x^2 + 256/625 over Rational Field
Elliptic Curve defined by y^2 = x^3 - x^2 - 8*x + 112 over Rational Field
Mapping from: CrvPln: C to CrvEll: E
with equations :
3/25*x^2 - 4/25*x*y + 1/25*y^2
-8/125*x^2 - 8/125*x*y + 8/125*y^2 + 8/25*x*z
-1/8*x^2 + 1/4*x*y - 1/8*y^2 + 1/4*x*z - 1/4*y*z
Elliptic curve isomorphism from: CrvEll: E to CrvEll: Em
Taking (x : y : 1) to (25/4*x + 2 : 125/8*y : 1)
1 true
Elliptic curve isomorphism from: CrvEll: Em to CrvEll: E
Taking (x : y : 1) to (4/25*x - 8/25 : 8/125*y : 1)
Mapping from: CrvEll: E to CrvPln: C
with equations :
25/2*$.1^3 - 35/2*$.1*$.2*$.3 + 25/4*$.2^2*$.3 + 56/5*$.1*$.3^2 - 8*$.2*$.3^2 +
    64/25*$.3^3
-5/2*$.1*$.2*$.3 + 75/4*$.2^2*$.3 + 8/5*$.1*$.3^2 - 24*$.2*$.3^2 + 192/25*$.3^3
$.1^2*$.3 - 15/2*$.1*$.2*$.3 + 24/5*$.1*$.3^2
and inverse
3/25*x^2 - 4/25*x*y + 1/25*y^2
-8/125*x^2 - 8/125*x*y + 8/125*y^2 + 8/25*x*z
-1/8*x^2 + 1/4*x*y - 1/8*y^2 + 1/4*x*z - 1/4*y*z

```
"""


def bremner_F1_point(k=1):
    coeffs_F1 = (Fraction(-1), Fraction(-30), Fraction(72))
    P = GeneralPoint(2, 4, coeffs_F1) * k
    X, Y = P.x, P.y

    x = (-21 + 4 * X) * (3 - 10 * X + 2 * X**2 + Y)
    y = 171 + 102 * X - 44 * X**2 - 69 * Y - 4 * X * Y
    z = (-21 + 4 * X) * (-15 + 5 * X + 2 * Y)

    A = (x**2 * y) - (x * y**2) + (3 * x * z**2) - (7 * y * z**2)
    B = z * (x**2 + 2 * x * y - y**2 + 10 * z**2)
    C = (x**2 * y) - (x * y**2) + (7 * x * z**2) - (3 * y * z**2)
    D = -z * (x**2 - 2 * x * y - y**2 - 10 * z**2)
    E = -z * ((x**2) + (y**2) + (2 * x * z) + (2 * y * z) - (10 * z**2))
    G = -z * ((x**2) + (y**2) - (2 * x * z) - (2 * y * z) - (10 * z**2))

    H_sq = (
        (x**4 * y**2)
        - (2 * x**3 * y**3)
        + (x**2 * y**4)
        + (10 * x**3 * y * z**2)
        - (20 * x**2 * y**2 * z**2)
        + (10 * x * y**3 * z**2)
        + (4 * x**3 * z**3)
        + (4 * x**2 * y * z**3)
        + (4 * x * y**2 * z**3)
        + (4 * y**3 * z**3)
        + (45 * x**2 * z**4)
        - (10 * x * y * z**4)
        + (45 * y**2 * z**4)
        - (40 * x * z**5)
        - (40 * y * z**5)
    )

    F_sq = (
        (x**4 * y**2)
        - (2 * x**3 * y**3)
        + (x**2 * y**4)
        + (10 * x**3 * y * z**2)
        - (20 * x**2 * y**2 * z**2)
        + (10 * x * y**3 * z**2)
        - (4 * x**3 * z**3)
        - (4 * x**2 * y * z**3)
        - (4 * x * y**2 * z**3)
        - (4 * y**3 * z**3)
        + (45 * x**2 * z**4)
        - (10 * x * y * z**4)
        + (45 * y**2 * z**4)
        + (40 * x * z**5)
        + (40 * y * z**5)
    )

    H = Fraction(math.isqrt(H_sq.numerator), math.isqrt(H_sq.denominator))
    F = Fraction(math.isqrt(F_sq.numerator), math.isqrt(F_sq.denominator))

    assert H**2 == H_sq
    assert F**2 == F_sq

    sol = [A, B, C, D, E, F, G, H]

    lcm = sp.lcm([t.denominator for t in sol])
    sol = [t * lcm for t in sol]
    gcd = sp.gcd(sol)
    sol = [int(t // gcd) for t in sol]

    return sol


def bremner_F2_point(k=5):
    coeffs = (Fraction(-1), Fraction(-8), Fraction(112))
    P = GeneralPoint(Fraction(12), Fraction(40), coeffs) * k

    X = Fraction(4, 25) * P.x - Fraction(8, 25)
    Y = Fraction(8, 125) * P.y

    x = (
        Fraction(25, 2) * X**3
        - Fraction(35, 2) * X * Y
        + Fraction(25, 4) * Y**2
        + Fraction(56, 5) * X
        - 8 * Y
        + Fraction(64, 25)
    )
    y = (
        Fraction(-5, 2) * X * Y
        + Fraction(75, 4) * Y**2
        + Fraction(8, 5) * X
        - 24 * Y
        + Fraction(192, 25)
    )
    z = X**2 - Fraction(15, 2) * X * Y + Fraction(24, 5) * X

    A = x**2 * y - x * y**2 + (3 * x - 7 * y) * z**2
    B = z * (x**2 + 2 * x * y - y**2 + 10 * z**2)
    C = x**2 * y - x * y**2 + (7 * x - 3 * y) * z**2
    D = -z * (x**2 - 2 * x * y - y**2 - 10 * z**2)
    E = -z * (x**2 + y**2 + 2 * (x + y) * z - 10 * z**2)
    F = (
        -2 * x**2 * y
        + 3 * x * y**2
        - y**3
        - (x**2 - 8 * x * y + 7 * y**2) * z
        - (x + 3 * y) * z**2
    )
    G = -z * (x**2 + y**2 - 2 * (x + y) * z - 10 * z**2)
    H = -(x**2) * y + x * y**2 + (6 * x * y - 2 * y**2) * z + (3 * x - 11 * y) * z**2

    sol = [A, B, C, D, E, F, G, H]

    lcm = sp.lcm([t.denominator for t in sol])
    sol = [t * lcm for t in sol]
    gcd = sp.gcd(sol)
    sol = [int(t // gcd) for t in sol]

    return sol


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


def is_quadruplet_E3(a, b, c, d):
    a_sq = a**2
    b_sq = b**2
    c_sq = c**2
    d_sq = d**2

    P = a_sq + b_sq
    if P % 2 == 1:
        return None, False
    if P != (c_sq + d_sq):
        return None, False
    P = P // 2
    Q = ((a_sq - b_sq) ** 2 + (c_sq - d_sq) ** 2) // 8
    R = ((a_sq - b_sq) ** 2 - (c_sq - d_sq) ** 2) ** 2 // 64
    # R = ((a_sq * b_sq) - (c_sq * d_sq)) ** 2 // 4
    return verify_skew_normal_squareish([P, Q, R])


class GaussianIntegersParameterization:
    def __init__(self, order):
        self.order = order
        self.perms = list(itertools.product([0, 1], repeat=order - 1))
        if self.order <= 1:
            return
        self.guide = list(itertools.product([0, 1], repeat=self.order))
        self.re_guide = [g for g in self.guide if sum(g) % 2 == 0]
        self.im_guide = [g for g in self.guide if sum(g) % 2 == 1]
        self.im_sign_matrix = (
            -(
                (
                    (
                        (np.array(self.im_guide)[:, 1:])
                        @ np.array(self.perms).transpose()
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
                        @ np.array(self.perms).transpose()
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
        if self.order > 1:
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
        else:
            rv = np.array(parameters[0]).reshape((2, 1))

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


def get_all_decomposable_primes_up_to(b, debug=True):
    ps = sp.primerange(b)
    return [2] + [p for p in tqdm.tqdm(ps, total=b, disable=not debug) if p % 4 == 1]


def get_all_decomposed_primes_up_to(b, debug=True):
    d_primes = get_all_decomposable_primes_up_to(b, debug=debug)
    return {p: decompose_prime(p) for p in tqdm.tqdm(d_primes, disable=not debug)}


def get_all_gaussian_integers_with_norm(N):
    rv = []
    for a in range(math.floor(math.sqrt(N)) + 1):
        for b in range(
            math.floor(math.sqrt(N - a**2)), min(a, int(math.sqrt(N - a**2))) + 1
        ):
            if a**2 + b**2 == N:
                rv.append((a, b))
    return rv


def factorize_gaussian_integer(a, b):
    norm = a**2 + b**2
    factors = sp.factorint(norm)
    integer_factors = {p: m for p, m in factors.items() if p % 4 == 3}
    integer_part = calculate_norm_from_factorization(integer_factors)
    decomposable_factors = {p: m for p, m in factors.items() if p % 4 != 3}
    num_primes_including_mult = sum(m for m in decomposable_factors.values())
    prime_decompositions = [
        decompose_prime(p) for p, m in decomposable_factors.items() for _ in range(m)
    ]
    gp = GaussianIntegersParameterization(num_primes_including_mult)
    logger.trace(
        f"\n {a=}, {b=}\n {norm=}\n {factors}\n {integer_factors=}\n {integer_part=}\n {decomposable_factors=}\n {prime_decompositions=}\n {num_primes_including_mult=}"
    )
    results = [
        tuple(t)
        for t in (
            math.isqrt(integer_part)
            * gp.get_pair_group_at(
                prime_decompositions, canonize_to_first_quadrant=True
            )
        )
        .transpose()
        .tolist()
    ]
    logger.trace(results)
    canonical_form = tuple(sorted((abs(a), abs(b))))
    i = results.index(canonical_form)
    conj_mask = (0,) + gp.perms[i]
    logger.trace(conj_mask)

    rec = math.prod(
        sp.ZZ_I(p[0], p[1] * (2 * do_conj - 1))
        for p, do_conj in zip(prime_decompositions, conj_mask)
    ) * math.isqrt(integer_part)
    assert canonical_form == tuple(sorted((abs(rec.x), abs(rec.y))))
    logger.trace(rec)
    return integer_part, prime_decompositions, conj_mask


def calculate_norm_from_factorization(factorization):
    rv = 1
    for p, m in factorization.items():
        rv *= p**m
    return rv


def find_quadruplets(points):
    points_a2b2 = [(p[0] ** 2) * (p[1] ** 2) for p in points]
    Q_m_dict = collections.defaultdict(list)
    for i, norm_i in enumerate(points_a2b2):
        for j, norm_j in enumerate(points_a2b2[:i]):
            Q_m_dict[norm_i + norm_j].append((i, j))
    return {
        k: ((points[v[0][0]], points[v[0][1]]), (points[v[1][0]], points[v[1][1]]))
        for k, v in Q_m_dict.items()
        if len(v) > 1
    }


def get_all_gaussian_integers_with_factored_norm(
    factors_list, return_metadata=False, prime_decompositions=dict()
):
    gaussian_primes = []
    j_p_ranges = []

    one_mult_p_representative = None
    if 2 in factors_list:
        assert factors_list[2] == 1
    for p, e_p in factors_list.items():
        gaussian_primes.append(p)
        if p not in prime_decompositions:
            prime_decompositions[p] = decompose_prime(p)
        if p == 2:
            j_p_ranges.append([1])
            continue
        assert p % 4 == 1, f"{p}"

        if (one_mult_p_representative is None) and e_p % 2 == 1:
            one_mult_p_representative = p
            j_p_ranges.append(list(range(1, e_p + 1, 2)))
        else:
            j_p_ranges.append(list(range(-e_p, e_p + 1, 2)))

    assert one_mult_p_representative is not None, factors_list
    rv = []
    j_ps_sequence = list(itertools.product(*j_p_ranges))

    for j_ps in j_ps_sequence:
        res_sp = 1
        for i, p in enumerate(gaussian_primes):
            d = prime_decompositions[p]
            d_sp = sp.ZZ_I(*d)
            d_conj_sp = sp.ZZ_I(*d[::-1])
            if j_ps[i] > 0:
                res_sp *= d_sp ** j_ps[i]
            else:
                res_sp *= d_conj_sp ** (-j_ps[i])
            res_sp *= p ** ((factors_list[p] - abs(j_ps[i])) // 2)
        rv.append(tuple(sorted([abs(res_sp.x), abs(res_sp.y)], reverse=True)))
    if return_metadata:
        return rv, j_ps_sequence, gaussian_primes, factors_list
    else:
        return rv


def validate_and_canonize_E4_sol(sol):
    sol = [abs(x) for x in sol]
    sol_pairs = [sorted(x, reverse=True) for x in zip(sol[::2], sol[1::2])]
    sol_quads = [sorted(x, reverse=True) for x in zip(sol_pairs[::2], sol_pairs[1::2])]
    sol_octs = [sorted(x, reverse=True) for x in zip(sol_quads[::2], sol_quads[1::2])]
    sol_quads = [x[i] for x in sol_octs for i in range(2)]
    sol_pairs = [x[i] for x in sol_quads for i in range(2)]
    sol = [x[i] for x in sol_pairs for i in range(2)]

    gcd = int(sp.gcd(sol))
    sol = [t // gcd for t in sol]
    radius = sol[0] ** 2 + sol[1] ** 2
    assert radius == sol[2] ** 2 + sol[3] ** 2
    assert radius == sol[4] ** 2 + sol[5] ** 2
    assert radius == sol[6] ** 2 + sol[7] ** 2
    Q_m = (sol[0] * sol[1]) ** 2 + (sol[2] * sol[3]) ** 2
    assert Q_m == (sol[4] * sol[5]) ** 2 + (sol[6] * sol[7]) ** 2
    if radius % 2 != 0:
        sol_pairs_mul2 = [(x[0] + x[1], x[0] - x[1]) for x in sol_pairs]
        return validate_and_canonize_E4_sol(
            [x[i] for x in sol_pairs_mul2 for i in range(2)]
        )

    return sol


def canonize_syndrome_pattern(p):
    x = np.array([list(_t) for _t in p.split(",")]) == "1"
    x = np.pad(x, [[0, 0], [1, 0]]).transpose()
    perms = [
        [0, 1, 2, 3],
        [1, 0, 2, 3],
        [0, 1, 3, 2],
        [1, 0, 3, 2],
        [2, 3, 0, 1],
        [2, 3, 1, 0],
        [3, 2, 0, 1],
        [3, 2, 1, 0],
    ]

    augs = []
    for perm in perms:
        x_p = x[perm]
        x_p ^= x_p[0:1, :]
        augs.append(
            tuple(
                np.sort(
                    np.sum(
                        x_p.astype(np.int32) * np.array([[8], [4], [2], [1]]), axis=0
                    )
                ).tolist()
            )
        )
    return str(sorted(augs)[0])


def analyze_E4_sol(sol, factorize=True, factorize_gaussian_integers=False):
    radius = sol[0] ** 2 + sol[1] ** 2
    Q_m = (sol[0] * sol[1]) ** 2 + (sol[2] * sol[3]) ** 2
    Q_a = sol[0] ** 4 + sol[1] ** 4 + sol[2] ** 4 + sol[3] ** 4
    L_2 = (sol[0] ** 2 - sol[1] ** 2) ** 2 + (sol[2] ** 2 - sol[3] ** 2) ** 2
    Q_b = L_2 - 4 * Q_m

    assert radius % 2 == 0
    L_1 = radius // 2
    assert L_2 % 8 == 0
    L_2 = L_2 // 8

    rv = dict()

    rv["A"] = sol[0]
    rv["B"] = sol[1]
    rv["C"] = sol[2]
    rv["D"] = sol[3]
    rv["E"] = sol[4]
    rv["F"] = sol[5]
    rv["G"] = sol[6]
    rv["H"] = sol[7]
    rv["lg2_coeff_size"] = len(bin(rv["A"])) - 2

    rv["L_1"] = L_1
    if factorize:
        rv["L_1_fact"] = sp.factorint(L_1)
    # rv["L_2"] = L_2
    # rv["L_2_fact"] = sp.factorint(L_2)
    # rv["Q_a"] = Q_a
    # rv["Q_a_fact"] = sp.factorint(Q_a)
    rv["Q_b"] = Q_b
    if factorize:
        rv["Q_b_fact"] = sp.factorint(Q_b)
    # rv["Q_m"] = Q_m
    # rv["Q_m_fact"] = sp.factorint(Q_m)

    if factorize_gaussian_integers:
        i_part_00, factors_00, conj_mask_00 = factorize_gaussian_integer(sol[0], sol[1])
        i_part_01, factors_01, conj_mask_01 = factorize_gaussian_integer(sol[2], sol[3])
        i_part_10, factors_10, conj_mask_10 = factorize_gaussian_integer(sol[4], sol[5])
        i_part_11, factors_11, conj_mask_11 = factorize_gaussian_integer(sol[6], sol[7])

        assert i_part_00 == i_part_01
        assert i_part_00 == i_part_10
        assert i_part_00 == i_part_11

        assert factors_00 == factors_01
        assert factors_00 == factors_10
        assert factors_00 == factors_11

        xor_mask = lambda x, y: tuple(xx ^ yy for xx, yy in zip(x, y))
        xorred_00 = xor_mask(conj_mask_00, conj_mask_00)
        xorred_01 = xor_mask(conj_mask_00, conj_mask_01)
        xorred_10 = xor_mask(conj_mask_00, conj_mask_10)
        xorred_11 = xor_mask(conj_mask_00, conj_mask_11)
        syndrome_pattern = sorted(
            list(
                set(
                    [
                        "".join([str(_x) for _x in x])
                        for x in zip(xorred_11, xorred_10, xorred_01)
                    ]
                )
            )
        )

        rv["syndrome_pattern"] = canonize_syndrome_pattern(",".join(syndrome_pattern))
    return rv
