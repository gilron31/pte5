# GPU-Accelerated Enumeration — Algorithm Draft

Working spec for the compiled enumerator behind `paper.md` §6. This is the
*real* enumerator; `code/native/` is a capability sandbox (GMP/CUDA/pybind
plumbing proven on toy kernels), not this.

Every formula below is verified in `drafts/verify_claims.py` or in the profiling
runs quoted in §2. **Nothing here is a size estimate pulled from the air** —
where a number is extrapolated, it says so.

---

## 0. What we are computing

Three searches, from `paper.md` §6:

| Search | Input | Output | §6 ref |
| --- | --- | --- | --- |
| **S-A** per-norm $E_4$ enumeration | admissible $L_1 \le B$ | the $E_4$ catalog | 6.1–6.2 |
| **S-B** top-collision | catalog + families | $E_5$, or exclusion to $B$ | 6.3 |
| **S-C** bottom-lift | catalog, $\mu \le 10^6$ | $E_5$, or null | 6.4 |
| **S-D** $\mathbb{F}_p$ tower sieve | $p \le 10^3$, depth $\le 5$ | densities, wheels, nonexistence branch | 6.6 |

S-A is the load-bearing one: Proposition 6.1 turns an exhaustive S-A catalog
into the headline exclusion. S-A is also where essentially all the compute goes.

---

## 1. The reduction we actually implement

Work at the **norm-$L_1$ level** ($L_1$ odd, canonical), not the root level.

**Points.** $U(L_1) = \{\, z = m+ni \in \mathbb{Z}[i] : N(z) = L_1 \,\}$ modulo
units and conjugation, canonized to the first octant $m \ge n \ge 0$.
$|U| = D = \lceil \tfrac12\prod_i (b_i+1) \rceil$ for $L_1 = \prod q_i^{b_i}$.

**The twist to roots.** $z \mapsto (1+i)z$ sends norm $L_1$ to norm $2L_1$; a
point $z = m+ni$ gives the root pair $(r,s) = (m{-}n,\, m{+}n)$ with
$r^2+s^2 = 2L_1$. The $(1{+}i)$ factor is applied once, at output time.

**The invariant.** Define
$$
\nu(z) := \mathrm{Re}(z^2) = m^2 - n^2 \;=\; r\cdot s .
$$
*Verified on the exemplar:* $z = 149794 + 103091i$, $N(z) = L_1$,
$\nu = 11810488155 = A\cdot B$. $\nu^2$ is invariant under the full unit group
and conjugation, so it is well defined on canonized points.

**The condition.** Four points $z_1,\dots,z_4 \in U$ split as $\{z_1,z_2\} \mid \{z_3,z_4\}$
form an $E_4$ iff
$$
\boxed{\;\nu_1^2 + \nu_2^2 \;=\; \nu_3^2 + \nu_4^2\;}
$$
(this is $L_{2,b}$; *verified*: $1207464164080651264050$ on both sides of the
exemplar) **and** the eight roots are distinct and nonzero (Proposition 5.4).

So **S-A is a disjoint-pair-sum collision problem on $D$ integers.** That is the
entire inner loop.

### 1.1 The product form (the algorithmic lever)

For a pair $\{z_1, z_2\}$ write $g = \gcd(z_1,z_2)$, $h = z_1/g$; then
$z_2 = \bar h g$ up to a unit, and
$$
L_{2,d}(\text{pair}) \;=\; 2\,\mathrm{Re}(g^4)\,\mathrm{Re}(h^4),
\qquad
L_{2,d} = 8L_1^2 - 8L_{2,b},
$$
both *verified* on the exemplar ($g = 165+467i$, $h = 259-450i$, product
$= -912832201535971887688 = L_{2,d}$). Since $L_1$ is fixed per norm, matching
$L_{2,b}$ ⟺ matching $L_{2,d}$, so we may use **either** key.

The consequence is structural: **the pair key is a product of two independent
factors**, one per side of a splitting of the prime multiset, and $\mathrm{Re}(\cdot^4)$
is invariant under units *and* conjugation. The pair keys are therefore not
"random sums" — they are products drawn from two precomputable tables. Whether
this yields a sub-quadratic collision search is **open** (see §5), but it is the
most promising lever for the champion regime, where the quadratic cost bites.

---

## 2. Measured workload shape

Exact sieve counts (not estimates) of admissible $L_1 \le B$ — products of split
primes only, with $D \ge 4$:

| $B$ | split-prime products | admissible ($D\ge4$) | pair-sums | avg pairs/norm | max $D$ |
| --- | --- | --- | --- | --- | --- |
| $10^6$ | 87,880 | 11,083 | 99,354 | 9.0 | 16 |
| $10^7$ | 814,181 | 131,437 | 1,339,509 | 10.2 | 24 |
| $10^8$ | 7,618,314 | 1,464,114 | 16,790,533 | 11.5 | 36 |

Growth is $\approx 11\times$ per decade in norms, $\approx 12.5\times$ in pairs.
**Extrapolated** to the campaign frontier:

| $B$ | admissible norms | pair-sums |
| --- | --- | --- |
| $10^{10}$ | $\sim1.7\times10^{8}$ | $\sim2.6\times10^{9}$ |
| $10^{12}$ | $\sim2\times10^{10}$ | $\sim4\times10^{11}$ |

### 2.1 The distribution is extremely top-heavy — this drives the design

At $B = 10^8$: $D=4$ is **75.5%** of admissible norms, and $D \in \{4,6,8\}$ is
**96.4%**. Two consequences:

- The dominant case is *tiny*. $D=4$ admits exactly **three** disjoint pair
  partitions ($12|34$, $13|24$, $14|23$) — three comparisons, no hash table, no
  shared memory, everything in registers.
- General $D$: the number of disjoint-pair-pairs is $D(D{-}1)(D{-}2)(D{-}3)/8$
  ($D{=}6\to45$, $D{=}8\to210$), versus $O(D^2)$ for hash-then-check-disjoint.
  Direct comparison wins up to $D \approx 8$–$10$, i.e. for ~97% of all norms.

**S-A is therefore not compute-bound; it is a throughput problem over ~$10^{10}$
tiny, highly uniform tasks.** Load balancing and norm *generation*, not FLOPs,
are the bottleneck.

# GILRO i think i might be able to suggest something here. If we find some way to schedule norms such that norms with a significat overlap in factors will be scheduled near each other. That way, some of the point generation might be amortized. And maybe norm generation can happen in the shared memory

---

## 3. Precision policy

$L_1 \le 10^{12} \Rightarrow m,n < 10^6$ (32-bit), $\nu \le L_1 < 2^{40}$
(64-bit), but $\nu_1^2+\nu_2^2 \le 2L_1^2 \approx 2\times10^{24}$ — **81 bits**.
Already at the exemplar ($L_1 \approx 3.3\times10^{10}$) the key is **71 bits**,
so 64-bit exact arithmetic is unusable above $L_1 \approx 2\times10^{9}$.

# GILRO why mersenne prime and not just $2^k$ for $k=64/32/16/8$?

**Policy: filter modulo the Mersenne prime $P = 2^{61}-1$, verify exactly on the
host.** Reduction is shift-and-add (no division); $\nu$ is exact in 64 bits, and
$\nu^2 \bmod P$ costs one 64×64→128 multiply plus a fold.

The property that matters for the exclusion claim:

> **Modular filtering has no false negatives.** Equal keys are always congruent,
> so a true $E_4$ can never be missed. Only false positives occur, at rate
> $2^{-61}$ per comparison — across the whole $B = 10^{12}$ campaign
> ($\sim4\times10^{11}$ comparisons) the expected count of spurious hits is
> $\sim2\times10^{-7}$. Every hit is re-verified in exact arithmetic anyway.

This is the right asymmetry for a null result (§6.8): the failure mode of a
filter is extra work, never a silent miss. A 32-bit filter would also be
sound but would produce $\sim10^{2}$ false positives at $10^{12}$ — still fine,
and worth having as the reduced-precision fast path.

---

## 4. Algorithms

### A1 — Factor base (CPU, one-time, reusable)

For every split prime $p \le P_{\max}$, store $(p, a, b)$ with $p = a^2+b^2$,
$a > b > 0$, via Tonelli–Shanks for $\sqrt{-1}$ then Cornacchia.

**Sizing.** An admissible norm needs $\ge 3$ split primes, the two smallest
being $\ge 5\cdot13 = 65$, so the largest prime factor can reach $B/65$. For
$B=10^{12}$ that is $P_{\max} \approx 1.5\times10^{10}$, i.e.
$\pi_1(P_{\max}) \approx 3\times10^{8}$ primes $\approx 4$ GB at 12 B/entry.
This table — not the arithmetic — is the real resource constraint. Options:
segment it (§A2), delta-encode $p$, or recompute per segment.

### A2 — Norm generation: **largest-prime-outer**

The single most important structural decision. Enumerate each admissible norm
*exactly once by its largest prime factor*:

```
for each split prime p (streamed in segments, ascending):
  for e = 1, 2, ... while p^e <= B:
    for each c in SmallTable with c <= B/p^e and maxprime(c) < p:
      emit norm n = c * p^e     with factorization and point set inherited
```

`SmallTable` = all split-prime products $\le B/p_{\min}$, precomputed once with
their point sets $U(c)$. For $p > 10^6$ the cofactor bound is $< 10^6$, so
`SmallTable` has only ~88k entries — small, shared, cache/constant-resident.

Why this ordering:

1. **Cornacchia is amortized.** Each prime's Gaussian decomposition is computed
   once and reused across all its multiples — $3\times10^8$ decompositions
   total instead of one per norm ($2\times10^{10}$).
2. **$U$ is inherited, not rebuilt.** $U(c\,p) = U(c)\mathfrak{p} \cup U(c)\bar{\mathfrak{p}}$:
   one Gaussian multiply per point instead of a full factor-by-factor rebuild.
3. **Perfect streaming shape.** The outer stream (primes) is sequential and
   segmentable; the inner table is small and reused — exactly the access
   pattern a GPU wants.
4. **No factoring, ever.** Norms arrive pre-factored by construction.

### A3 — Point-set construction

```
U(n) = dedup_first_octant( U(c) * 𝔭^j * conj(𝔭)^(e-j),  j = 0..e )
```

Multiplicative doubling with the conjugate trick (4 int-multiplies yield both
$z\mathfrak{p}$ and $z\bar{\mathfrak{p}}$):
$x{=}m a,\ y{=}m b,\ z{=}n a,\ w{=}n b \Rightarrow (x{-}w)+(z{+}y)i$ and
$(x{+}w)+(z{-}y)i$. Branching over *each copy* of a repeated prime is what
realizes exponents $0 \le c \le b_i$ — the non-squarefree coverage that
Proposition 5.1 requires and that every exclusion claim rests on.

### A4 — The meet (hot loop), by $D$-bucket

Bucket norms by $D$ and launch a kernel per bucket, so every lane does identical
work with zero divergence.

- **$D = 4$** (~75%): compute $\kappa_i = \nu_i^2 \bmod P$; test the 3
  partitions. Registers only.
- **$D = 6, 8$** (~21%): 45 / 210 disjoint partitions, unrolled, registers only.
- **$D \ge 10$** (~3%): compute $D(D{-}1)/2$ pair keys, sort or hash in shared
  memory, report equal keys, then check disjointness.
- **Champion regime** ($D \sim 2^{19}$, §6.2(ii)): $D^2/2 \approx 1.4\times10^{11}$
  keys cannot be materialized. Generate pair keys in tiles and keep only those
  whose hash falls in a residue window of width $1/R$; repeat $R$ times. At
  $R = 10^3$ each pass retains ~$1.4\times10^8$ keys (~1 GB), and the total is
  $R \cdot D^2/2 \approx 1.4\times10^{14}$ key evaluations — hours, not days.

### A4′ — Product-form meet (**open**, §1.1)

Use $\kappa_{\text{pair}} = \mathrm{Re}(g^4)\mathrm{Re}(h^4)$ and search for
*multiplicative* collisions over splittings of the prime multiset. If this
admits a divisor-based or meet-in-the-middle attack it changes the champion
regime from quadratic to something better. Unproven; prototype before relying on
it.

### A5 — Exact verification and tower recovery (CPU, rare)

Every filter hit: recompute the key in exact 128-bit; recover
$(L_1,L_2,L_3,\sqrt{L_4})$ by the polynomial peeling of §3; assert the degree-16
identity $\prod(x^2-r_k^2) = (((x^2{-}L_1)^2{-}L_2)^2{-}L_3)^2 - L_4$ and the
§4 ladder divisibilities. Cheap (hits are rare) and it is the correctness spine.

### A6 — Canonization and dedup

Scale out common factors; rotate to the first octant; normalize to $L_1$ odd;
sort pairs descending. Emits the catalog key and the syndrome pattern —
**pending the canonization lemma** (§5 TODO), which for repeated split primes
must record a *multiset* over buckets. Until that lemma lands, the syndrome
census is provisional; the catalog itself is not affected.

### A7 — Top-collision join (S-B)

Hash-join the catalog on $(L_1,L_2,L_3)$; any collision with distinct $L_4$ is
an $E_5$. This is an external sort / hash-join over $\sim10^{8}$–$10^{10}$
records, i.e. a database problem, not an arithmetic one. Radix-partition by a
hash of the key so each partition fits in RAM.

### A8 — Bottom-lift join (S-C)

For each catalogued $E_4$ with roots $r'_0..r'_7$, seek $(\mu, L)$ with
$L \pm \mu r'_k$ both squares for all $k$. Per $k$ this is a divisor enumeration
of $2\mu r'_k$; join the eight on the shared $L$. Parity wheel: under the odd-$L_1$
normalization the roots are odd, forcing $\mu \equiv 0 \pmod 4$, which kills
$\mu = 1,2,3$ immediately.

### A9 — $\mathbb{F}_p$ tower sieve (S-D)

Do **not** brute-force $p^5$. Use the §3 recursion as a DP over levels: a
depth-$d$ tower over $\mathbb{F}_p$ is a depth-$(d{-}1)$ tower in $u$-space plus
the congruum condition, so each level is an $O(p^2)$ pass over residue sets.
Outputs: local densities $\sigma_p$ for §8, sieve wheels for A2, and the
jackpot branch — a single $p$ with zero depth-5 towers is an unconditional
nonexistence theorem.

---

## 5. Open algorithmic questions

1. **Is A4′ sub-quadratic?** The product structure of the pair key is real and
   verified; exploiting it is not yet designed. Highest-value open item.
2. **Factor-base memory.** 4 GB at $B=10^{12}$ vs 8 GB of device memory on the
   current card. Segment, compress, or compute Cornacchia on-GPU per segment?
3. **Where does the $D$-bucketing boundary sit** between the register-only
   comparison kernel and the hash kernel? Measure; the crossover is predicted
   around $D \approx 8$–$10$ but depends on occupancy.
4. **Is the $(1+i)$ twist worth deferring?** Working at norm $L_1$ halves the
   operand size; confirm no solution is lost at the boundary.
5. **Champion-norm selection.** Which norms maximize $D$ per unit of work — pure
   primorial-style products of the smallest split primes, or something tuned to
   the ladder congruences of §4?

---

## 6. Correctness obligations (from §6.8)

Non-negotiable for an exclusion claim, and cheapest to build in from the start:

- Two independent implementations (the Python enumerator in `code/old/` is the
  natural oracle for the compiled one).
- **Family self-test**: re-find every member of Bremner's two families,
  *including non-squarefree $L_1(k)$* — the test that exposes an A3 coverage hole.
- Synthetic injection: plant known $E_4$s in a segment, confirm recovery.
- Per-norm counts cross-checked against the A9 model; deviations are bugs, not
  discoveries.
- Pre-registered $B$, checkpointed logs, resumable segments.

# GILRO general remarks:

1. Be clear about where does norm generation and point generation happens within the gpu (or were you postponing it to the implementation design?)
2. How are we  addressing multiplicities and dedups? Essentially for every norm with a square split prime (actually every multiplicity larger  that 2), a third of the time it will just be a scaled up version of the norm divided by that square. Are we bypassing these situations somehow?
