# GPU-Accelerated Enumeration — Algorithm Draft (rev. 2)

Working spec for the compiled enumerator behind `paper.md` §6. This is the
*real* enumerator; `code/native/` is a capability sandbox (GMP/CUDA/pybind
plumbing proven on toy kernels), not this.

Every formula is verified in `drafts/verify_claims.py` or in the experiments
quoted inline. **No number here is a guess** — where a figure is extrapolated,
it says so.

> **§5 is superseded by `drafts/gpu_implementation_plan.md` (rev. 3).** §1–§4 —
> the reduction, the injectivity lemma, the workload shape, the algorithms — are
> unchanged and remain the reference. What changed is the *strategy*: the CPU
> campaign was built and run, S-A to the design bound turned out to be minutes,
> and the ranking of what the GPU is for is now measured rather than guessed.
> Three numbers in §4 also moved: A4's champion-regime cost is $1.4\times10^{11}$
> key evaluations, not $1.4\times10^{14}$ (window by key *value*, not by hash);
> A9's tower sieve is $O(p^4/64)$ word-ops per prime with an explicit recursion;
> and A9's "one prime with zero depth-5 towers is an unconditional nonexistence
> theorem" **is false as stated** — it holds only for splitting *with
> multiplicity*, and the distinct-root reading gives zero at 25 of the 29 primes
> below 113 while a known $E_4$ fails that same test at 16 of them.

---

> **Reading guide.** §1 reduces the whole $E_4$ search to a single one-line
> condition, and §1.3 shows that condition is a *product* equality — the
> structural fact everything else exploits. §2 measures how much work there
> actually is and discovers the shape is lopsided in a way that dictates the
> kernel design. §3–§4 are the engineering consequences. If you read only two
> things, read §1.2 (the injectivity lemma, which deletes an entire step) and
> §2.2 (the $D=4$ dominance, which decides the kernels).

---

## 0. What we are computing

| Search | Input | Output | §6 ref |
| --- | --- | --- | --- |
| **S-A** per-norm $E_4$ enumeration | admissible $L_1 \le B$ | the $E_4$ catalog | 6.1–6.2 |
| **S-B** top-collision | catalog + families | $E_5$, or exclusion to $B$ | 6.3 |
| **S-C** bottom-lift | catalog, $\mu \le 10^6$ | $E_5$, or null | 6.4 |
| **S-D** $\mathbb{F}_p$ tower sieve | $p \le 10^3$, depth $\le 5$ | densities, wheels, nonexistence | 6.6 |

S-A is load-bearing: Proposition 6.1 turns an exhaustive S-A catalog into the
headline exclusion. S-A is also where essentially all the compute goes.

---

## 1. The reduction we actually implement

We work at the **norm-$L_1$ level** ($L_1$ odd, canonical), never at the root
level. The reason is size: roots live at norm $2L_1$, so staying downstairs
keeps every operand one bit smaller and lets the $(1{+}i)$ twist be applied once
at output.

**Points.** $U(L_1) = \{\, z = m+ni \in \mathbb{Z}[i] : N(z) = L_1 \,\}$ modulo
units and conjugation, canonized to the first octant $m \ge n \ge 0$. For
$L_1 = \prod_i q_i^{b_i}$,
$$
D := |U(L_1)| = \Big\lceil \tfrac12 \textstyle\prod_i (b_i+1) \Big\rceil .
$$

**The twist to roots.** $z \mapsto (1+i)z$ maps norm $L_1$ to norm $2L_1$; the
point $z = m+ni$ becomes the root pair $(r,s) = (m{-}n,\, m{+}n)$ with
$r^2+s^2 = 2L_1$.

**The invariant.** Define
$$
\nu(z) := \mathrm{Re}(z^2) = m^2 - n^2 \;=\; r\cdot s .
$$
*Verified on the exemplar:* $z = 149794 + 103091i$ has $N(z) = L_1$ and
$\nu = 11810488155 = A\cdot B$. Because $\mathrm{Re}((iz)^2) = -\nu$ and
$\mathrm{Re}(\bar z^2) = \nu$, the square $\nu^2$ is invariant under the full
unit group *and* conjugation — so it is well defined on canonized points.

**The condition.** Four points split as $\{z_1,z_2\} \mid \{z_3,z_4\}$ form an
$E_4$ iff
$$
\boxed{\;\nu_1^2 + \nu_2^2 \;=\; \nu_3^2 + \nu_4^2\;}
$$
(*verified*: both sides $=1207464164080651264050$ on the exemplar), together
with the nondegeneracy of Proposition 5.4.

**So S-A is a disjoint-pair-sum collision on $D$ integers.** That is the entire
inner loop.

### 1.2 Injectivity of $\nu$ — this deletes a whole step

**Lemma.** *$\nu$ is injective on canonized points of a fixed norm.*

*Proof.* $m^2+n^2 = L_1$ and $m^2-n^2 = \nu$ give $m^2 = (L_1+\nu)/2$,
$n^2 = (L_1-\nu)/2$, so $(m,n)$ is determined by $\nu$. $\square$

*Checked on 300 random norms; no collision ever.* Two consequences, both
practical:

1. **Disjointness is automatic.** If two *distinct* pairs share a key, they
   cannot share a point: pairs $\{a,b\}$ and $\{a,c\}$ with equal key force
   $\nu_b^2 = \nu_c^2$, hence $b = c$. So every key collision is already a
   genuine four-point configuration — **the kernel needs no disjointness
   test at all.** (Nondegeneracy of the *roots* is still Proposition 5.4, but
   that is a cheap output-side check.)
2. $\nu$ is a perfect hash of the point, so a point can be stored as its $\nu$
   alone — one 64-bit word instead of two coordinates.

### 1.3 The product form (the algorithmic lever)

Take a pair $\{z_1,z_2\}$. Let $g$ be the product of the Gaussian primes where
the two points *agree*, $h$ where they *differ*, so $z_1 = gh$ and $z_2 = g\bar h$
up to a unit. Writing $\mathrm{Re}(z^4) = 2\nu^2 - L_1^2$ and
$z_1^4 + z_2^4 = g^4(h^4 + \bar h^4)$:

$$
\boxed{\;\nu_1^2+\nu_2^2 \;=\; \mathrm{Re}(g^4)\,\mathrm{Re}(h^4) \;+\; L_1^2\;}
$$

*Verified on 12,088 pairs across random norms: zero mismatches.* Since $L_1$ is
fixed per norm, the additive constant cancels and:

> **The $E_4$ condition is exactly a product equality**
> $\;\mathrm{Re}(g_1^4)\mathrm{Re}(h_1^4) = \mathrm{Re}(g_2^4)\mathrm{Re}(h_2^4)$
> between two different splittings of the prime multiset.

This is `paper.md` §5's "$L_{2,d}$ prism" in computational form. It matters
because the key is no longer an opaque sum: it is a product of two
independently-tabulable factors, each invariant under units and conjugation, and
each itself a product of two $\mathbb{Z}[\sqrt2]$ norm-form values
(Theorem 5.3: $\mathrm{Re}(g^4) = (m^2{-}n^2{-}2mn)(m^2{-}n^2{+}2mn)$). See §4A4′
for what this does and does not buy.

---

## 2. Measured workload shape

Exact sieve counts (not estimates) of admissible $L_1 \le B$ — products of split
primes only, with $D \ge 4$:

| $B$ | split-prime products | admissible ($D\ge4$) | pair-sums | avg pairs/norm | max $D$ |
| --- | --- | --- | --- | --- | --- |
| $10^6$ | 87,880 | 11,083 | 99,354 | 9.0 | 16 |
| $10^7$ | 814,181 | 131,437 | 1,339,509 | 10.2 | 24 |
| $10^8$ | 7,618,314 | 1,464,114 | 16,790,533 | 11.5 | 36 |

Growth $\approx11\times$/decade in norms, $\approx12.5\times$ in pairs.
**Extrapolated:**

| $B$ | admissible norms | pair-sums |
| --- | --- | --- |
| $10^{10}$ | $\sim1.7\times10^{8}$ | $\sim2.6\times10^{9}$ |
| $10^{12}$ | $\sim2\times10^{10}$ | $\sim4\times10^{11}$ |

At $B = 10^{12}$ the richest admissible norm is
$5^2\!\cdot\!13\!\cdot\!17\!\cdot\!29\!\cdot\!37\!\cdot\!41\!\cdot\!53\!\cdot\!61 \approx 7.9\times10^{11}$,
giving $D_{\max} = 192$. **$D$ stays small.** That is the single most important
engineering fact in this document.

### 2.1 The distribution is lopsided — this decides the kernels

At $B = 10^8$: $D=4$ is **75.5%** of admissible norms; $D \in \{4,6,8\}$ is
**96.4%**.

- $D = 4$ admits exactly **three** disjoint pair-partitions ($12|34$, $13|24$,
  $14|23$). Three comparisons. No hash table, no shared memory, no
  disjointness check (§1.2) — pure registers.
- For general $D$ the number of disjoint pair-pairs is
  $D(D{-}1)(D{-}2)(D{-}3)/8$ ($D{=}6\to45$, $D{=}8\to210$) versus $O(D^2)$ for
  hash-and-collide. Direct comparison wins to $D \approx 8$–$10$, i.e. for ~97%
  of all norms.

**S-A is therefore a throughput problem over $\sim10^{10}$ tiny, uniform tasks,
not a compute-bound one.** Generation and scheduling — not arithmetic — set the
pace.

### 2.2 Scheduling by factor overlap (this is what makes it GPU-shaped)

The natural amortization is the **cofactor DFS tree**: order norms so that
consecutive ones differ by a single prime. Then a child's point set comes from
its parent by *one* Gaussian multiply per point,
$U(cq) = U(c)\mathfrak{q} \cup U(c)\bar{\mathfrak{q}}$, instead of a rebuild
from scratch.

Assign **one block per subtree**. The ancestor path's point sets live in shared
memory; a block walks its subtree emitting thousands of norms with *no global
memory traffic at all* except the (essentially nonexistent) hits.

Sizing: depth $\le 8$ and $D \le 192$ at $B=10^{12}$, so the whole path costs
$8 \times 192 \times 8\,\text{B} = 12$ KB — comfortable against 48–96 KB of
shared memory, and far smaller for the typical subtree.

This is what resolves the arithmetic-intensity problem. A naive "CPU generates,
GPU meets" split moves 32 bytes to do ~15 integer ops — **0.5 ops/byte, roughly
80× PCIe-bound, and slower than not using the GPU at all.** Generating inside
the block instead raises intensity by three orders of magnitude and makes the
device-resident design (§5, option B) the right one.

---

## 3. Precision policy

$L_1 \le 10^{12} \Rightarrow m,n < 10^6$ (32-bit) and $\nu \le L_1 < 2^{40}$
(64-bit), but the key $\nu_1^2+\nu_2^2 \le 2L_1^2$ needs **81 bits** — already
**71 bits** at the exemplar. So exact 64-bit fails above $L_1 \approx 2\times10^9$
and we must filter.

**Why not simply truncate mod $2^k$?** Because the keys carry forced 2-adic
structure, so the low bits are nearly information-free. Since $m^2+n^2 = L_1$ is
odd, $\nu$ is odd and $\nu^2 \equiv 1 \pmod 8$; more sharply, Dilcher's Remark 3
says the $L_2$-difference is always divisible by $2^7\cdot3^2$. Measured
directly on 4,000 random key differences at the norm-$L_1$ level:

$$
\min v_2 = 6, \qquad \min v_3 = 2
\qquad\Longrightarrow\qquad 2^6\cdot3^2 = 576 \ \big|\ (\text{every difference}).
$$

A modulus of $2^k$ therefore throws away its bottom 6 bits:

| modulus | effective bits | false positives over $4\times10^{11}$ comparisons |
| --- | --- | --- |
| $2^{16}$ | 10 | $4\times10^{8}$ — useless |
| $2^{32}$ | 26 | $6\times10^{3}$ — usable as a first-stage filter |
| $2^{64}$ | 58 | $10^{-6}$ — **free and ample** |

**Policy (revised): use natural `uint64` overflow as the primary filter.** It
costs nothing, and 58 effective bits is far more than enough. A prime modulus
would recover the lost 6 bits but buys nothing we need; drop the earlier
Mersenne proposal. Keep $2^{32}$ in reserve as a cheap pre-filter if register
pressure ever demands it.

The property that matters for the exclusion claim is unchanged:

> **Modular filtering has no false negatives.** Equal keys are always congruent,
> so a true $E_4$ can never be missed. The failure mode is extra verification
> work, never a silent miss — exactly the asymmetry a null result needs (§6.8).

---

## 4. Algorithms

### A1 — Factor base

For each split prime $p \le P_{\max}$, obtain $p = a^2+b^2$ via Tonelli–Shanks
for $\sqrt{-1}$ then Cornacchia.

An admissible norm needs $\ge3$ split primes and the two smallest are $\ge 5\cdot13 = 65$,
so $P_{\max} = B/65 \approx 1.5\times10^{10}$ at $B=10^{12}$ — about
$3\times10^8$ primes, 4 GB stored. **Do not store it.** Cornacchia is ~150
integer ops, so the whole base is $\sim4.5\times10^{10}$ ops — seconds of GPU
time. Recompute per segment; the memory problem disappears.

### A2 — Norm generation: largest-prime-outer

Enumerate each admissible norm *exactly once, by its largest prime factor*:

```
for each split prime p (segmented, ascending):          # device-resident
  for e = 1, 2, ... while p^e <= B:
    for each c in SmallTable with c <= B/p^e and maxprime(c) < p:
      emit n = c * p^e     # factorization and point set inherited
```

`SmallTable` = split-prime products $\le B/p$, held as a DFS tree (§2.2). For
$p > 10^6$ the cofactor bound is $<10^6$, so it has only ~88k nodes — small
enough to be resident.

Why this ordering:

1. **Cornacchia is amortized** — once per prime ($3\times10^8$) instead of once
   per norm ($2\times10^{10}$).
2. **$U$ is inherited, not rebuilt** — one Gaussian multiply per point.
3. **No factoring, ever** — norms arrive pre-factored by construction.
4. **It streams**: primes are sequential and segmentable, the tree is reused.

### A3 — Where the work physically lives

Answering the "which level of the GPU" question explicitly:

| Stage | Location | Notes |
| --- | --- | --- |
| Prime segment sieve + Cornacchia | device, one kernel per segment | output stays in global memory |
| `SmallTable` DFS tree | device global, read-mostly | ~88k nodes, built once |
| Subtree walk (norm generation) | **per block**, state in **shared memory** | one block owns one subtree |
| Point-set $U$ construction | **shared memory**, incremental from parent | 1 Gaussian mult/point |
| $\nu$ values, keys, the meet | **registers** | $D \le 192$; $D\le8$ is 97% of cases |
| Hits | global, atomically appended | rare; verified on host |

Gaussian multiply uses the conjugate trick — four 32-bit multiplies yield *both*
$z\mathfrak{q}$ and $z\bar{\mathfrak{q}}$:
$x{=}ma,\ y{=}mb,\ z{=}na,\ w{=}nb \Rightarrow (x{-}w)+(z{+}y)i$ and $(x{+}w)+(z{-}y)i$.

Branching over *each copy* of a repeated prime realizes exponents
$0 \le c \le b_i$ — the non-squarefree coverage Proposition 5.1 requires, and on
which every exclusion claim rests.

### A3.1 — Multiplicities, imprimitive points, and duplicates

This needs care, and it is easy to get silently wrong.

For $L_1 = p^{\,b}M$ with $p \nmid M$, a point with $p$-exponent $c$ (i.e.
$\mathfrak{p}^c\bar{\mathfrak{p}}^{\,b-c}$) is divisible by $p$ **iff**
$1 \le c \le b-1$. So the imprimitive fraction of $U$ is $(b{-}1)/(b{+}1)$:

| $b$ | 1 | 2 | 3 | 4 |
| --- | --- | --- | --- | --- |
| imprimitive share of $U$ | 0 | **1/3** | 1/2 | 3/5 |

and that imprimitive family is *exactly* $p\cdot U(M)$ — a scaled copy of a
smaller norm's point set. *Verified* at $L_1 = 5^2\cdot13\cdot17$: $|U| = 6$,
imprimitive $=2 = |U(13\cdot17)|$, and the two sets coincide after scaling.

The consequence must be stated precisely, because the obvious rule is wrong:

- A quadruple with **all four** points imprimitive at $p$ is $p^2 \times$ a
  solution at $M$ — a **duplicate**, already found when $M$ was processed.
  Rejecting these is required, or every solution reappears at $p^2M$, $q^2M$, …
  and the catalog counts (and any density fit) are corrupted.
- A quadruple with **some but not all** points imprimitive is **primitive and
  must be kept**. This is precisely the S3 case of `critical_review.md`: one
  pair uses the rational factor $\mathfrak{p}\bar{\mathfrak{p}} = p$, another
  uses $\mathfrak{p}^2$. Discarding it would reopen the coverage hole that
  Proposition 5.1 exists to close.

**Rule:** reject a hit iff $\gcd(z_1,z_2,z_3,z_4) \ne 1$ — implemented as a
bitmask test "all four have $1 \le c_p \le b_p-1$" for some $p$, using the
exponent vectors already carried by A2. Cost: a few ALU ops, on hits only.

How much does this cost us in wasted work? Pairs lying entirely inside the
imprimitive family re-derive the norm-$M$ search; for $b=2$ that is
$|U(M)|^2 / |U(p^2M)|^2 \approx 1/9$ of the pair work at such norms. And
non-squarefree norms are a minority of norms but carry a disproportionate share
of pair work (larger $D$):

| $B$ | non-squarefree share of norms | of pair work |
| --- | --- | --- |
| $10^6$ | 21.5% | 37.2% |
| $10^7$ | 19.0% | 36.2% |
| $10^8$ | 17.3% | 35.4% |

So the recoverable saving is ~4% overall — not worth optimizing. **Handle this
for correctness, not for speed.**

### A4 — The meet (hot loop), by $D$-bucket

Bucket by $D$ so every lane does identical work. Since $D(cq) = 2D(c)$, the
bucket is known *before* generation — sort `SmallTable` by $D$ once and the
batches are homogeneous for free, with no stream compaction.

- **$D = 4$** (~75%): compute $\kappa_i = \nu_i^2$ (mod $2^{64}$); test 3
  partitions. Registers only.
- **$D = 6,8$** (~21%): 45 / 210 partitions, unrolled, registers only.
- **$D \ge 10$** (~3%): compute $D(D{-}1)/2$ keys, sort or hash in shared
  memory, report equal keys. **No disjointness check** (§1.2).
- **Champion regime** ($D \sim 2^{19}$): $D^2/2 \approx 1.4\times10^{11}$ keys
  cannot be materialized. Generate in tiles, keep only keys whose hash lands in
  a residue window of width $1/R$, repeat $R$ times. At $R=10^3$ each pass
  retains $\sim1.4\times10^8$ keys ($\sim$1 GB) for a total of
  $\sim1.4\times10^{14}$ key evaluations — hours, not days.

### A4′ — Product-form meet: what the exploration found

Using §1.3, a pair key is $\mathrm{Re}(g^4)\mathrm{Re}(h^4)$ over a splitting of
the prime multiset into agree-part and differ-part. Findings:

**Confirmed.** The identity holds exactly (12,088 pairs, 0 mismatches), and
every $\mathrm{Re}(g^4)$ factors into two $\mathbb{Z}[\sqrt2]$ norm-form values
whose odd primes are all $\equiv\pm1 \pmod 8$ (Theorem 5.3, re-verified here).

**What it buys.**

1. *Arithmetic:* for a fixed differ-set $S$ the keys form a rank-1 grid
   $A_i \cdot B_j$ over precomputed tables $R(S^c,\cdot)$ and $R(S,\cdot)$, so a
   key costs **one multiply and two loads** instead of two squarings and an add
   — a constant-factor win of perhaps $2\times$, no more.
2. *A structural filter:* a collision is an equality of **prime multisets**.
   Theorem 5.3 confines those primes to $\equiv\pm1\pmod 8$, halving the
   candidate prime set, and each $\mathrm{Re}(g^4)$ is a product of two values
   of size $\sim L_1^{1/2}$ — small enough to factor cheaply. Bucketing pairs by
   $v_\pi$ for a few small admissible $\pi$ is then a *structural*, not random,
   pre-filter.

**What it does not buy (yet).** No asymptotic gain. Summing over all differ-sets
$S$ still costs $\sum_S 2^{s-2} = 2^{2s-3}$ — identical to the naive pair count.
The grid structure converts the search into "find $A_i/A_k = B_l/B_j$", which is
still quadratic unless the ratio structure is exploited.

**The live research direction.** Because collisions are equalities of
factorizations over a *restricted* prime set, the natural framing is
index-calculus-like: factor the $R$-values over a basis of admissible primes,
represent each as an exponent vector, and search for 4-term relations by linear
algebra rather than by enumeration. Whether the relation density makes this
competitive is open, and is the single highest-value algorithmic question here.

**Verdict:** A4′ is a *champion-regime* tool. In the exhaustive regime $D \le 8$
and the naive key is already three ALU ops — there is nothing to optimize.

### A5 — Exact verification and tower recovery (host, rare)

On every filter hit: recompute in exact 128-bit; recover
$(L_1,L_2,L_3,\sqrt{L_4})$ by the peeling of §3; assert the degree-16 identity
$\prod(x^2-r_k^2) = (((x^2{-}L_1)^2{-}L_2)^2{-}L_3)^2 - L_4$, the §4 ladder
divisibilities, and primitivity (A3.1). This is the correctness spine and costs
nothing, since hits are almost nonexistent.

### A6 — Canonization and dedup

Scale out common factors; rotate to the first octant; normalize $L_1$ odd; sort
pairs descending. Emits the catalog key and the syndrome pattern — **pending the
canonization lemma** (`paper.md` §5 TODO), which for repeated split primes must
record a *multiset* over buckets. Until it lands the syndrome census is
provisional; the catalog itself is unaffected.

### A7 — Top-collision join (S-B)

Hash-join the catalog on $(L_1,L_2,L_3)$; any collision with distinct $L_4$ is
an $E_5$. A database problem, not an arithmetic one: radix-partition by a hash
of the key so each partition fits in RAM.

### A8 — Bottom-lift join (S-C)

For each catalogued $E_4$ with roots $r'_0..r'_7$, seek $(\mu,L)$ with
$L \pm \mu r'_k$ both squares for all $k$; per $k$ a divisor enumeration of
$2\mu r'_k$, joined on the shared $L$. Parity wheel: odd roots force
$\mu \equiv 0 \pmod 4$, killing $\mu = 1,2,3$ outright.

### A9 — $\mathbb{F}_p$ tower sieve (S-D)

Do **not** brute-force $p^5$. Use the §3 recursion as a DP over levels: a
depth-$d$ tower over $\mathbb{F}_p$ is a depth-$(d{-}1)$ tower in $u$-space plus
the congruum condition, so each level is an $O(p^2)$ pass. Outputs: local
densities $\sigma_p$ for §8, sieve wheels for A2, and the jackpot branch — one
prime with zero depth-5 towers is an unconditional nonexistence theorem.

---

## 5. Implementation strategy

The arithmetic-intensity analysis of §2.2 rules out the naive hybrid. Two
viable designs remain:

- **(A) Split by regime.** CPU (C++/OpenMP) owns the exhaustive campaign; GPU
  owns champion norms and the $\mathbb{F}_p$ sieve, where intensity is enormous.
  Simplest; $B = 10^{10}$–$10^{11}$ is comfortable on CPU alone.
- **(B) Device-resident.** Norm generation, point construction and the meet all
  live on the GPU; nothing crosses PCIe but the prime stream and hits. The
  subtree-in-shared-memory scheme of §2.2 is what makes this work.

Recommendation: **build A first** — it is the natural oracle for B (§6) and gets
a real catalog quickly — then port to B if the campaign targets $10^{12}$.

Hardware on hand: RTX 4060 Laptop, 8 GB, sm\_89, CUDA 12.9. `double` runs at
1/64 rate, so every kernel is integer-only.

## 6. Open questions

1. **Is A4′'s index-calculus framing competitive?** The only route to breaking
   the champion regime's quadratic wall. Highest value, genuinely uncertain.
2. **Where does the $D$-bucket crossover sit** between register-only comparison
   and shared-memory hashing? Predicted $D \approx 8$–$10$; occupancy-dependent,
   so measure.
3. **Subtree granularity.** How deep a subtree per block balances shared-memory
   footprint against scheduling overhead?
4. **Champion-norm selection.** Maximize $D$ per unit work — primorial-style
   products of the smallest split primes, or tuned to the §4 ladder congruences?

## 7. Correctness obligations (from §6.8)

- Two independent implementations — `code/old/enumeration_v1.py` is the natural
  oracle for the compiled one.
- **Family self-test**: re-find every member of Bremner's two families,
  *including non-squarefree $L_1(k)$* — the test that exposes an A3/A3.1 hole.
- Synthetic injection: plant known $E_4$s in a segment; confirm recovery.
- Per-norm counts cross-checked against the A9 model; deviations are bugs, not
  discoveries.
- Pre-registered $B$, checkpointed logs, resumable segments.

---

## Revision notes (rev. 1 → rev. 2)

| Comment | Resolution |
| --- | --- |
| *Schedule norms with overlapping factors together; generate in shared memory* | Adopted as the core design — §2.2 (cofactor DFS tree, one block per subtree) and §A3. It is also what fixes the arithmetic-intensity problem and makes the device-resident option viable. |
| *Why a Mersenne prime and not $2^k$?* | **You were right.** §3 rewritten: measured $\min v_2 = 6$ on key differences, so $2^k$ loses only 6 bits — `uint64` overflow gives 58 effective bits, free. Mersenne proposal dropped. |
| *Where do norm/point generation happen on the GPU?* | Was postponed; now explicit in §A3 (table of stage → memory level). |
| *How are multiplicities and dedups handled? Isn't a third of a $p^2$-norm just a scaled copy?* | **Correct, and rev. 1 glossed over it.** New §A3.1: imprimitive share is exactly $(b{-}1)/(b{+}1)$ (= 1/3 at $b=2$), the imprimitive family *is* $p\cdot U(M)$, and the rule is subtle — reject only when **all four** points are imprimitive; keep partially-imprimitive quadruples, which are exactly the S3 primitive case. Quantified: ~4% recoverable work, so it is a correctness fix, not an optimization. |
| — | New §1.2: $\nu$ is injective, so key collisions are automatically disjoint — the disjointness check is removed from every kernel. |
