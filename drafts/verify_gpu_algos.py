"""Reproduces every measured number in drafts/gpu_algorithms.md (rev. 2).

Sections:
  1. nu injectivity  -> disjointness is automatic (doc 1.2)
  2. A4-prime product form  key = Re(g^4)Re(h^4) + L1^2  (doc 1.3, A4-prime)
  3. Re(g^4) prime content, Thm 5.3 filter strength
  4. modulus choice: min v_2 / v_3 of key differences (doc 3)
  5. multiplicities / imprimitive points / duplicates (doc A3.1)
  6. admissible-norm and D-distribution profile (doc 2)
"""

import itertools, math, random
from sympy import factorint, primerange


def gmul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def gconj(z):
    return (z[0], -z[1])


def gnorm(z):
    return z[0] * z[0] + z[1] * z[1]


def canon(z):
    a, b = abs(z[0]), abs(z[1])
    return (max(a, b), min(a, b))


def split_prime(p):
    for a in range(1, math.isqrt(p) + 1):
        b2 = p - a * a
        b = math.isqrt(b2)
        if b > 0 and b * b == b2:
            return (max(a, b), min(a, b))
    raise ValueError(p)


def nu(z):
    return z[0] * z[0] - z[1] * z[1]  # Re(z^2)


def Re4(z):  # Re(z^4) = (m^2-n^2)^2 - 4m^2n^2
    m, n = z
    return (m * m - n * n) ** 2 - 4 * m * m * n * n


def build_U(primes):
    """primes = list of split primes (with multiplicity). Return canonized points + sign vectors."""
    pts = {}
    qs = [split_prime(p) for p in primes]
    for signs in itertools.product([0, 1], repeat=len(qs)):
        z = (1, 0)
        for q, s in zip(qs, signs):
            z = gmul(z, q if s == 0 else gconj(q))
        pts.setdefault(canon(z), signs)
    return pts


print("=" * 72)
print(
    "1. nu is INJECTIVE on canonized points  =>  key collision is automatically disjoint"
)
print("=" * 72)
ok = True
tested = 0
sp = [p for p in primerange(3, 200) if p % 4 == 1]
for trial in range(300):
    k = random.randint(3, 6)
    ps = random.sample(sp, k)
    U = build_U(ps)
    nus = [nu(z) for z in U]
    tested += 1
    if len(set(nus)) != len(nus):
        ok = False
        print("  COLLISION at", ps)
print(f"  {tested} random norms, nu injective on U every time: {ok}")
print("  (proof: m^2+n^2=L1 and m^2-n^2=nu determine m,n uniquely)")

print()
print("=" * 72)
print("2. A4' PRODUCT FORM at the norm-L1 level")
print("=" * 72)
print(
    "  claim: key_b(pair) = nu1^2+nu2^2 = Re(g^4)*Re(h^4) + L1^2,  g=agree-part, h=differ-part"
)
bad = 0
checked = 0
for trial in range(200):
    k = random.randint(3, 6)
    ps = random.sample(sp, k)
    qs = [split_prime(p) for p in ps]
    L1 = 1
    for p in ps:
        L1 *= p
    for e1 in itertools.product([0, 1], repeat=k):
        for e2 in itertools.product([0, 1], repeat=k):
            if e1 >= e2:
                continue
            g = (1, 0)
            h = (1, 0)
            for q, a, b in zip(qs, e1, e2):
                if a == b:
                    g = gmul(g, q if a == 0 else gconj(q))
                else:
                    h = gmul(h, q if a == 0 else gconj(q))
            z1 = gmul(g, h)
            z2 = gmul(g, gconj(h))
            if canon(z1) == canon(z2):
                continue
            key_b = nu(canon(z1)) ** 2 + nu(canon(z2)) ** 2
            checked += 1
            if key_b != Re4(g) * Re4(h) + L1 * L1:
                bad += 1
    if trial > 20:
        break
print(f"  checked {checked} pairs, mismatches: {bad}")
print(
    "  => E4 condition is EXACTLY  Re(g1^4)Re(h1^4) == Re(g2^4)Re(h2^4)   (a PRODUCT equality)"
)

print()
print("=" * 72)
print("3. Re(g^4) prime content (Thm 5.3 mechanism) -- how rigid is the filter?")
print("=" * 72)
vals = []
for trial in range(60):
    k = random.randint(2, 4)
    ps = random.sample(sp, k)
    qs = [split_prime(p) for p in ps]
    g = (1, 0)
    for q in qs:
        g = gmul(g, q)
    v = abs(Re4(g))
    if v > 1:
        vals.append(v)
allp = set()
bad8 = 0
for v in vals:
    f = factorint(v)
    for p in f:
        if p != 2:
            allp.add(p)
            if p % 8 not in (1, 7):
                bad8 += 1
print(f"  {len(vals)} Re(g^4) values; odd prime factors all == +-1 mod 8: {bad8==0}")
print(
    f"  distinct odd primes seen: {len(allp)}; density of admissible primes below max:"
)
mx = max(allp)
adm = [p for p in primerange(3, min(mx, 20000)) if p % 8 in (1, 7)]
print(
    f"    primes<20000 with p==+-1 mod 8: {len(adm)} of {len(list(primerange(3,20000)))}  (~1/2 of odd primes)"
)
print("  => Thm 5.3 halves the candidate prime set; a real but not decisive filter.")

print()
print("=" * 72)
print("4. GILRO Q: why not modulus 2^k?  (Dilcher Rmk 3: 1152 | key difference)")
print("=" * 72)
v2s = []
v3s = []
for trial in range(4000):
    k = random.randint(3, 5)
    ps = random.sample(sp, k)
    U = list(build_U(ps))
    if len(U) < 4:
        continue
    a, b, c, d = random.sample(U, 4)
    diff = (nu(a) ** 2 + nu(b) ** 2) - (nu(c) ** 2 + nu(d) ** 2)
    if diff == 0:
        continue
    n = abs(diff)
    e2 = 0
    while n % 2 == 0:
        n //= 2
        e2 += 1
    e3 = 0
    while n % 3 == 0:
        n //= 3
        e3 += 1
    v2s.append(e2)
    v3s.append(e3)
print(f"  {len(v2s)} random key-differences at the norm-L1 level")
print(f"  min v_2 = {min(v2s)}   (mean {sum(v2s)/len(v2s):.1f})")
print(f"  min v_3 = {min(v3s)}   (mean {sum(v3s)/len(v3s):.1f})")
k2 = min(v2s)
k3 = min(v3s)
print(f"  => every difference is divisible by 2^{k2} * 3^{k3} = {2**k2 * 3**k3}")
print(f"  A modulus of 2^k therefore wastes its bottom {k2} bits:")
for K in (16, 32, 64):
    print(
        f"     mod 2^{K}: effective discrimination ~2^{K-k2}  "
        f"(FP rate {2.0**-(K-k2):.2e}; over 4e11 cmps -> {4e11*2.0**-(K-k2):.3g} false positives)"
    )
print("  A prime modulus coprime to 6 loses nothing.")


import itertools, math, numpy as np
from sympy import primerange


def gmul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def gconj(z):
    return (z[0], -z[1])


def canon(z):
    a, b = abs(z[0]), abs(z[1])
    return (max(a, b), min(a, b))


def split_prime(p):
    for a in range(1, math.isqrt(p) + 1):
        b2 = p - a * a
        b = math.isqrt(b2)
        if b > 0 and b * b == b2:
            return (max(a, b), min(a, b))
    raise ValueError(p)


print("=" * 72)
print("A. Structure of U(p^b * M): which points are IMPRIMITIVE (divisible by p)?")
print("=" * 72)
print("  point with p-exponent c (i.e. q^c qbar^(b-c)) is divisible by p iff 1<=c<=b-1")
for b in (1, 2, 3, 4):
    tot = b + 1
    imp = max(0, b - 1)
    print(
        f"   b={b}: c in 0..{b}, imprimitive c-values={imp}/{tot} = {imp/tot:.3f}"
        + ("   <-- GILRO's 'a third'" if b == 2 else "")
    )

# concrete: L1 = 5^2 * 13 * 17
ps = [5, 5, 13, 17]
qs = [split_prime(p) for p in ps]
U = {}
for signs in itertools.product([0, 1], repeat=len(qs)):
    z = (1, 0)
    for q, s in zip(qs, signs):
        z = gmul(z, q if s == 0 else gconj(q))
    U.setdefault(canon(z), signs)
L1 = 5 * 5 * 13 * 17
Ulist = list(U)
imp = [z for z in Ulist if z[0] % 5 == 0 and z[1] % 5 == 0]
print(
    f"\n  L1 = 5^2*13*17 = {L1}:  |U| = {len(Ulist)}, imprimitive (5|z) = {len(imp)} -> {len(imp)/len(Ulist):.3f}"
)
# the imprimitive family is exactly 5 * U(13*17)
U2 = {}
for signs in itertools.product([0, 1], repeat=2):
    z = (1, 0)
    for q, s in zip([split_prime(13), split_prime(17)], signs):
        z = gmul(z, q if s == 0 else gconj(q))
    U2.setdefault(canon(z), signs)
scaled = {canon((5 * z[0], 5 * z[1])) for z in U2}
print(
    f"  imprimitive family == 5*U(13*17): {set(imp)==scaled}   (|U(13*17)|={len(U2)})"
)

nu = lambda z: z[0] * z[0] - z[1] * z[1]
P = len(Ulist) * (len(Ulist) - 1) // 2
Pi = len(imp) * (len(imp) - 1) // 2
print(
    f"\n  pairs total={P}, pairs entirely inside imprimitive family={Pi} ({Pi/P:.1%})"
)
print(
    "  A quadruple with BOTH pairs imprimitive is p^2 x (a solution at M) -- a DUPLICATE."
)
print(
    "  A quadruple with only ONE pair imprimitive is PRIMITIVE and must be kept (the S3 case)."
)

print()
print("=" * 72)
print("B. How much of the campaign is non-squarefree? (exact sieve)")
print("=" * 72)


def profile(B):
    acc = np.ones(B + 1, dtype=np.int64)
    bad = np.zeros(B + 1, dtype=bool)
    sqfree = np.ones(B + 1, dtype=bool)
    s = np.ones(B + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(B**0.5) + 1):
        if s[i]:
            s[i * i :: i] = False
    for p in np.flatnonzero(s):
        p = int(p)
        if p == 2 or p % 4 == 3:
            bad[p::p] = True
        else:
            pe, e = p, 1
            while pe <= B:
                sl = acc[pe::pe]
                acc[pe::pe] = sl // e * (e + 1)
                pe *= p
                e += 1
        if p * p <= B:
            sqfree[p * p :: p * p] = False
    good = ~bad
    good[0] = good[1] = False
    adm = good & (acc >= 8)
    n_adm = int(adm.sum())
    n_nsf = int((adm & ~sqfree).sum())
    D = (acc + 1) // 2
    pairs_all = int((D[adm] * (D[adm] - 1) // 2).sum())
    pairs_nsf = int((D[adm & ~sqfree] * (D[adm & ~sqfree] - 1) // 2).sum())
    return n_adm, n_nsf, pairs_all, pairs_nsf


for e in (6, 7, 8):
    B = 10**e
    a, n, pa, pn = profile(B)
    print(
        f"  B=1e{e}: admissible={a:,}  non-squarefree={n:,} ({n/a:.1%})   "
        f"pair-work in nsf norms: {pn/pa:.1%}"
    )
print("\n  => non-squarefree norms are a MINORITY of norms but carry a")
print("     disproportionate share of pair work (larger D).")


def profile(B):
    # acc[n] = prod over split primes p|n of (v_p(n)+1); bad[n] = has factor 2 or p=3 mod 4
    acc = np.ones(B + 1, dtype=np.int64)
    bad = np.zeros(B + 1, dtype=bool)
    sieve = np.ones(B + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(B**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    primes = np.flatnonzero(sieve)
    for p in primes:
        p = int(p)
        if p == 2 or p % 4 == 3:
            bad[p::p] = True
        else:  # split prime: fold in (e+1)/e per prime power
            pe, e = p, 1
            while pe <= B:
                s = acc[pe::pe]
                acc[pe::pe] = s // e * (e + 1)
                pe *= p
                e += 1
    good = ~bad
    good[0] = good[1] = False
    adm = good & (acc >= 8)  # D>=4 needs prod(b+1)>=8
    nadm = int(adm.sum())
    D = (acc[adm] + 1) // 2
    pairs = int((D * (D - 1) // 2).sum())
    return int(good.sum()) - 1, nadm, pairs, int(D.max()), np.bincount(D)[:12]


for e in (6, 7, 8):
    B = 10**e
    tot, adm, pairs, dmax, hist = profile(B)
    print(
        f"B=1e{e}: split-products={tot:,}  admissible(D>=4)={adm:,}  pair-sums={pairs:,}  maxD={dmax}"
    )
    print(f"    avg pairs/norm={pairs/adm:.1f}   D-hist[4..11]={list(hist[4:12])}")
    sys.stdout.flush()
