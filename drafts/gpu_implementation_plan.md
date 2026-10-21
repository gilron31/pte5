# GPU Implementation Plan (rev. 3)

**What this supersedes.** `gpu_algorithms.md` §5 ("Implementation strategy") and
`implementation_design.md` §4 ("GPU critical path") both organise GPU work around
S-A: a K1/K2/K3 kernel pipeline with milestones M1–M4 ending in a device-resident
enumerator. That ordering was set before any campaign existed and is now
obsolete. **The mathematics of `gpu_algorithms.md` §1–§4 is unchanged** and
remains the reference; what changes is where the device goes and in what order.

**What this document is.** An answer to one question — *given a finished,
measured, oracle-grade CPU enumerator, what is the GPU actually for?* — followed
by sized designs for the two searches where the answer is "something".

> **Update.** The cost models below have since been replaced by measurements:
> see **`drafts/performance_ceiling.md`** and the `code/native/bench/` suite.
> The rankings held, the S-E numbers came in *better* than projected (138 s, not
> "minutes"), and one uncosted risk in §3.3 turned out to be free. Rows marked
> **[measured]** carry a real figure; everything else is still a model.

Hardware of record: RTX 4060 Laptop (sm_89, 8 GB, PCIe 4.0 ×8), CUDA 12.9,
against an AMD Ryzen 7 8845HS (8 cores / 16 threads, AVX-512).

---

## 0. The measurement that invalidated the old plan

The old plan's premise was `implementation_design.md` §3: a projected 16.5
core-hours for the direct meet at $B = 10^{12}$, which made S-A look like the
compute problem and the GPU like its answer. Bucketing took the projection to
~1 core-hour; then the campaign was built and run.

| | measured |
| --- | --- |
| $B = 10^{10}$, complete, 6 threads | **5.68 s** enumeration |
| admissible norms | 165,632,402 |
| segments | 4,323,837 |
| factor base (mmap of the $10^{12}$ artifact) | 0.23 s |
| catalogued, all certified | 8 |
| known-solution audit | 8/8 recovered |
| segment re-verification | 500 segments, 12,376 norms, 0 disagreements |

Extrapolating norms at the measured $11.1\times$/decade gives $\sim2\times10^{10}$
at the design bound, so at the observed rate S-A to $B = 10^{12}$ is **~12
minutes on six threads**, or 25–40 minutes under the 2–3× cache pessimism §3 of
the design doc already budgets. (Wall clock on this machine swings up to 2× from
frequency scaling — see `implementation_design.md` §9.4 — so read these as an
order of magnitude, not a stopwatch.)

> **S-A to the design bound is a coffee break on the CPU already sitting here.**
> Every milestone in the old GPU plan targeted a problem that no longer exists.

There is a second consequence, less obvious and more useful. The old plan had the
GPU and CPU cross-checking each other as co-equals. That is no longer the
relationship: **the CPU campaign is the oracle**, with an independent naive
`meet_oracle`, an audit module carrying known solutions from Bremner's families
and `code/old`, seeded segment re-verification, and a byte-identical catalog at
$10^{10}$ that agrees with the old catalogue in both directions. Any device code
is now checked *against* that, which makes a GPU port much cheaper to trust — and
means a port is no longer justified by "we need a second implementation".

---

## 1. Where the GPU earns its keep — ranked, with costs

| Search | CPU cost with today's code | GPU value | Verdict |
| --- | --- | --- | --- |
| **S-A** to $10^{12}$ | 12–40 min, 6 threads | ~10×, saves half an hour | **No** |
| **S-A** to $10^{14}$ | $\sim2\times10^{12}$ norms; 1–3 days, 16 threads | hours instead of days | **Only if $B$ moves** |
| **S-B** top-collision join | catalog-sized hash join; the $10^{10}$ catalog is 8 rows | — | **No** |
| **S-C** bottom-lift | unsized — never measured | unknown | **Measure on CPU first** |
| **S-D** $\mathbb{F}_p$ tower sieve | $O(p^4/64)$ word-ops/prime (§2); feasible to $p\approx10^3$ | **[measured]** rotate runs at 1.1e11 word-op/s at $p{=}65536$ but only **1.5e10 at $p{=}1024$** — the device is at its worst exactly where S-D lives | **Second**, CPU first — now on two independent grounds |
| **S-E** champion meet | $1.4\times10^{11}$ pair keys (§3, corrected down from $1.4\times10^{14}$) | **[measured]** collision detection **138 s** by hash insert, 163 s by radix sort; exact 256-bit windowing **1.4 s** | **First, and the clearest win** |

Two entries in that table moved by orders of magnitude while writing it, both
downwards, and both because an algorithm improved rather than because a machine
did. That is the pattern this project keeps repeating — the factor base, the
meet's bucketing, and now S-E — and it is the reason the ordering below puts
*measurement and algorithm work* ahead of every kernel.

---

## 2. S-D — the tower sieve, sized and corrected

> **[resolved] The model question this section opens is settled, and the answer
> changed the ring.** It is not $\mathbb{F}_p$ but $\mathbb{Z}/p^k$ with a cap on
> node valuations; the two models §2.2 weighs are the two extreme caps at $k=1$.
> `drafts/local_model_plan.md` is the argument and the measurements,
> `paper.md` §6.6 and Proposition 4.3 are the write-up, and
> `code/native/src/cpu/tower/` is the implementation. Read §2.1–§2.3 below as
> the sizing exercise it was; the recursion survives verbatim, the ring does not.

`paper.md` §6.6 asks, for each prime $p$ and depth $3 \le d \le 5$, for the count
of $(L_1,\dots,L_d) \in \mathbb{F}_p^d$ with $E_d$ splitting completely over
$\mathbb{F}_p$. It wants three products: local densities $\sigma_p$ for the
singular series of §8, sieve wheels for S-A, and **the nonexistence branch — a
single prime with zero depth-5 towers would prove unconditionally that no $E_5$
exists over $\mathbb{Z}$.** The third of those is false as stated (§2.2), and its
corrected form is a statement about $k \to \infty$ that no finite run settles.

### 2.1 The recursion, and why it is not $p^5$

$E_n(x) = E_{n-1}(x)^2 - L_n$ with $E_0(x) = x$. Define

$$T_d = \{\, c \in \mathbb{F}_p : E_d(x) - c \text{ splits into linear factors} \,\},
\qquad A_d = \{\, t \in T_d : -t \in T_d \,\}.$$

Since $E_d(x) - c = E_{d-1}(x)^2 - (L_d + c)$, a root needs $L_d + c$ to be a
square $t^2$, and then $E_d - c = (E_{d-1}-t)(E_{d-1}+t)$ splits iff both factors
do. Hence

$$\boxed{\,T_d = \{\, t^2 - L_d : t \in A_{d-1} \,\},\qquad T_0 = \mathbb{F}_p\,}$$

and $E_d$ splits iff $0 \in T_d$, i.e. **iff $L_d$ lies in the image of squaring
on $A_{d-1}$**. So for a fixed prefix $(L_1,\dots,L_{d-1})$ the number of $L_d$
completing it is $|\{t^2 : t \in A_{d-1}\}|$ — read off in $O(p)$ rather than by
retesting each $L_d$.

*Verified against brute-force polynomial factorization over $\mathbb{F}_p$ for
$p = 3,5,7,11,13$ and $d = 1,2,3$ — every prefix, exact set equality, plus the
counting corollary.* Script: `drafts/verify_gpu_plan.py`, check 1.

**The state collapses, and that is the whole algorithm.** The DP state is the set
$T$, not the tuple, and distinct reachable sets are few:

| $p$ | distinct states at level 1 / 2 / 3 / 4 |
| --- | --- |
| 13 | 13 / 26 / 14 / 1 |
| 31 | 31 / 931 / 497 / 32 |
| 71 | 71 / 4,971 / 7,527 / 72 |
| 113 | 113 / 3,277 / 13,674 / 1,809 |

States grow like $\sim p^2$, not $p^4$. With $T$ held as a $p$-bit mask, one
transition is $O(p/64)$ words, so a level costs $|states| \times p \times p/64
\approx p^4/64$ word-ops and $p^3$ bits of state. At $p = 10^3$ that is
$1.5\times10^{10}$ word-ops and 125 MB — **a CPU job of minutes.** At $p = 10^4$
the state is 125 GB and the approach dies; the ceiling is $p \sim$ a few
$\times10^3$ on this hardware regardless of who computes it.

### 2.2 The landmine: the obvious local model proves a false theorem

Running the DP with the natural nondegeneracy condition — require $2^d$
**distinct** roots, which is what an honest tower has, implemented as
$T^*_d = \{t^2 - L_d : t \in A^*_{d-1},\, t \ne 0\}$ (the $t \ne 0$ is exactly
disjointness of the $\pm t$ fibers) — gives:

| $p$ | 3 | 5 | 7 | … | 73 | **79** | 83 | **89** | **97** | 101 | 103 | 107 | 109 | **113** |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| nondegenerate depth-5 towers | 0 | 0 | 0 | 0 | 0 | 78 | 0 | 11 | 207 | 0 | 0 | 0 | 0 | 371 |

Twenty-five of the twenty-nine primes below 113 admit **zero** nondegenerate
depth-5 towers. Read naively that is the jackpot: an unconditional proof that no
$E_5$ exists over $\mathbb{Z}$, from $p = 3$.

**It is not, and the inference is wrong.** An integer $E_5$ has $32$ distinct
integer roots, but distinct integers may be *congruent mod $p$*. The reduction is
therefore a tower that splits **with multiplicity**, not one with distinct roots,
so $\sigma_p^{\text{distinct}} = 0$ excludes nothing. Checked on the §7 exemplar
$E_4$, whose 16 roots are distinct over $\mathbb{Z}$:

> the 16 roots collapse to fewer residues at **16 of the 29 primes $\le 113$** —
> only 3 distinct residues mod 5, 4 mod 7, 8 mod 13, 4 mod 37.

A real solution, in hand, fails the distinct-root local test at more than half of
these primes. Any nonexistence certificate must use the multiplicity-allowing
model — and that model has no zeros ($p=3$: 3 towers, $p=13$: 46, $p=23$: 199).

The correction cuts the other way too. Read the same table as a ratio: at the 25
primes where the distinct count is 0, the multiplicity count is not — $p=13$ has
46, $p=73$ has 5,122, $p=109$ has 15,715 — so **every tower those primes count
has a repeated root.** The mechanism is visible in the recursion: once $A$
collapses to $\{0\}$ the only completion is $L_d = 0$, i.e. $E_d = E_{d-1}^2$.
So at small $p$ the multiplicity-allowing $\sigma_p$ is a census of degeneracies,
and is no better as a density than the distinct model was as an obstruction.

> **Neither obvious local model is right.** Distinct-roots is too strong (proves
> a false theorem); with-multiplicity is too weak (counts only degeneracies at
> small $p$). Fixing this is a mathematics question, not an engineering one, and
> it must be answered before any $\sigma_p$ enters the singular series of §8 or
> any wheel prunes S-A. Likely shape of the answer: count mod $p^k$ for $k$ large
> enough that root collisions are resolved, or track the $p$-adic valuations of
> pairwise root differences. **This is now the top open question in the plan**,
> and it is CPU work — the DP that answers it runs in seconds at these $p$.

### 2.3 What the kernel looks like, once the model is settled

Genuinely GPU-shaped, and the shape is unusual enough to state:

- **State**: $p$-bit masks, $p/64$ words. One block owns one state; threads own
  words. All of $T$, $A$, and the squaring image are bitwise.
- **Transition**: for each of $p$ values of $L_d$, $T' = \{t^2 - L_d : t \in A\}$
  — a *permuted scatter* of the set $\{t^2 : t \in A\}$ by a rotation of index.
  Precompute $S = \{t^2 : t \in A\}$ once per state ($O(p)$), then every $L_d$ is
  a cyclic shift of $S$ by $L_d$: $p$ transitions for $O(p^2/64)$ words, and the
  shift is `funnelshift` — one instruction.
- **Dedup**: hash the $p$-bit mask (one 64-bit fingerprint per block, then verify
  on collision), insert into a device hash table with `atomicCAS`. This is the
  only irregular part and the only place the GPU version is harder than the CPU
  one.
- **Output**: per-level tuple counts, accumulated with one `atomicAdd` per state.

Occupancy is bounded by the mask: $p/8$ bytes per resident state. At $p = 10^3$
that is 125 B — trivial. At $p = 10^4$, 1.25 KB, still fine; it is host memory
for the state *table* that binds, not shared memory.

**Verdict: build it on the CPU first.** The DP is small, the model question is
unsettled, and a GPU kernel that computes the wrong quantity fast is worse than
no kernel. The GPU version's only purpose is to push $p$ past $10^3$, and that is
worth doing only if $\sigma_p$ turns out to converge slowly.

> **[updated] The CPU build happened and it moved the device verdict.** Two
> changes. **(1) The modulus is $p^k$, not $p$.** §5 ranks S-D behind S-E partly
> because P5 measured the bitset rotate at 1.5e10 word-op/s at $p = 1024$ against
> 1.1e11 at 65536 — "the device is at its worst exactly where S-D lives". The
> corrected model lives at $p^k \approx 2\times10^3$ to $8\times10^4$, which is
> the regime where P5 measured the device at its **best**. **(2) The dedup is
> the wrong thing to accelerate.** Measured: at $3^8$, level 2 performs
> $6561 \times 6561 = 4.3\times10^7$ interns that produce only $32{,}806$
> distinct states — **99.92% hits**, because the $m$ rotations of one state are
> re-discovering one translation orbit $m$ times. The fix is algorithmic
> (`local_model_plan.md` §5.1: canonicalise once per state instead of interning
> $m$ times, a factor $m/64$), and it deletes the `atomicCAS` table that §2.3
> calls "the only irregular part". A device port should follow that rewrite, not
> precede it.

---

## 3. S-E — the champion meet, with a 1000× algorithmic correction

The champion regime picks single hand-chosen norms with enormous $D$ and asks
for collisions inside them; it supports *discovery*, not exclusion. $D \sim 2^{19}$
needs ~20 distinct split primes, so $L_1 \approx 2.5\times10^{36}$ — far outside
any exhaustive bound, which is the point.

### 3.1 The correction

Both drafts size this as: $D^2/2 \approx 1.37\times10^{11}$ pair keys cannot be
materialised (1.1 TB), so keep only keys whose *hash* lands in a residue window
of width $1/R$ and repeat $R$ times — at $R = 10^3$, $\sim1.4\times10^{14}$ key
evaluations.

That recomputes every key $R$ times because a hash window is unpredictable. Window
by the key's **value** instead. With $sq[\,]$ sorted ascending, the pairs whose
sum falls in $[K_r, K_{r+1})$ are, for each $i$, a *contiguous* run of $j$ —
found by two binary searches. Each pair then belongs to exactly one window:

$$\text{total work} \;=\; \binom{D}{2} \;+\; R\,D\log D
\;=\; 1.37\times10^{11} + 1\times10^{10},$$

against $R\binom{D}{2} = 1.37\times10^{14}$. *Verified* on a random 400-element
corpus with $R=7$: 79,800 pairs enumerated, all distinct, exactly
$\binom{400}{2}$ — no pair lost, none repeated (`drafts/verify_gpu_plan.py`, check 4).

**Three orders of magnitude, from replacing a hash window with a sorted-order
window.** It applies to the CPU implementation identically.

### 3.2 What that leaves, and why it is still the best GPU candidate

Work per window: enumerate $\sim1.4\times10^8$ pairs, detect equal keys among
them. Key generation is now trivial (one add from an L2-resident table); the cost
is the collision detection — a sort or a hash over ~1.4×10⁸ items, i.e. a
*bandwidth* problem, which is exactly what a GPU is for. Per window at 16 B per
entry (64-bit truncated key + packed $(i,j)$, both indices $< 2^{19}$): 2.2 GB.
Take $R = 4\times10^3$ for ~550 MB per window, comfortably inside 8 GB.

**[measured]** — `bench/gpu/bench_gpu_primitives` P3 and P4, on the 4060:

| route | rate | S-E's $1.374\times10^{11}$ pairs |
| --- | --- | --- |
| CUB radix sort, 8+8 B pairs | 8.4e8 key/s (54 GB/s) | **163 s** |
| `atomicCAS` open addressing | 1.0e9 insert/s | **138 s**, including the table clear |

**Build the hash, not the sort**: it is faster and needs a quarter of the memory,
so the window count $R$ can be chosen for convenience rather than to fit 2 GB of
temp storage. The CPU equivalent pays a cache miss per key on a table far larger
than L3. Note the device clock *sagged to 1.09–1.37 GHz* during the large sorts
(power-limited laptop), so these are conservative.

### 3.3 The gap nobody has costed: the champion regime is outside the precision policy

`precision.h` and `gpu_algorithms.md` §3 are both sized for $L_1 \le 10^{12}$:
coordinates 20 bits, $\nu$ 40, the exact key 81. At $L_1 \approx 2.5\times10^{36}$
the coordinates are 61 bits, $\nu$ is 122, and $sq = \nu^2$ is **~243 bits — four
64-bit limbs.** The `uint64` filter still filters soundly (equal keys stay
congruent, so still no false negatives; at 58 effective bits and $1.4\times10^8$
keys per window the birthday load is ~0.03 false positives per window, ~130 over
the run — free to reject on the host).

**[measured] The multi-limb cost, flagged here as an uncosted risk, is nil.**
P6 puts a fixed 4-limb add plus lexicographic compare at 1.02e11 op-pairs/s, so
S-E's $1.374\times10^{11}$ exact windowing comparisons cost **1.4 s** against the
138 s of collision detection. Use exact 4-limb comparisons and stop worrying
about it. The remaining requirements stand:

- the **windowing comparison must be exact**, or a pair lands in the wrong window
  and is silently lost — the one failure mode Part IV of the design doc says is
  the only one that matters. Four-limb compare, not `double`;
- $sq[\,]$ is $D \times 32$ B = 16.8 MB, which is L2/global, not shared;
- point-set *construction* at this size is multi-limb Gaussian arithmetic, not
  the `int32` path.

`src/gpu/bigint/` and `src/cpu/bigint/` already exist in the sandbox with a
shared limb convention (32-bit limbs, least significant first) — they were built
as capability demos and this is the first search that actually needs them. What
is missing is a fixed-width 4-limb type with comparison and addition, which is
much less than the general `kMaxLimbs = 32` schoolbook kernel already there.

### 3.4 What it did when it was built — **[measured]**

G3 and G4 are built (`code/native/src/core/`, `src/cpu/champion/`,
`src/gpu/champion/`, `apps/champion_demo.cc`). The sized champion has been run:

| | |
| --- | --- |
| champion norm | $L_1 = 2386442377984630986454048784862092465$ (the 20 smallest split primes) |
| $L_1$ / $\nu$ / key | 121 / 121 / **243 bits** — as §3.3 projected |
| $D$ | 524,288 |
| pairs enumerated | **137,438,691,328** = $\binom{2^{19}}{2}$ exactly |
| point-set construction | 1.0 s (GMP; the `int32` path stops at 20-bit coordinates) |
| device time | **280 s**, 4.91e8 pair/s, 5,794 windows |
| candidates | 4,836 — **all** rejected by the exact 256-bit re-check |
| collisions | **none** |

Three things that measurement settled:

1. **The false-positive model is right.** The device compares a signature of
   table_log + 23 bits, not the key; at 2.37e7 pairs per window into 2^26 slots
   the model says ~0.5 spurious candidates per window and the run produced 0.83.
   Rejecting them costs one 256-bit compare each, so the filter/arbiter split of
   `precision.h` carries over to the champion regime unchanged.
2. **The run is exactly half host, and the kernel is exactly at its ceiling.**
   Instrumented (`ChampionStats::plan_seconds` / `device_seconds`):

   | | |
   | --- | --- |
   | device | **138.8 s** → 9.9e8 pair/s in-kernel |
   | host window planning | **139.3 s** → 50% of the run, device idle throughout |

   The device half reproduces P4's 1.0e9 insert/s to two digits, which is the
   cost model closing: the plan projected 138 s from a primitive benchmark and
   the search took 138.8 s. So the kernel is *done* until the other half is
   fixed — window planning is an $O(D)$ sweep per window, 5,794 of them over
   524,288 points, and the device idles through all of it. Cheapest fix is not
   a device-side planner but a second stream: plan window $N{+}1$ while window
   $N$ runs, double-buffering the range arrays, which is ~2× for a few dozen
   lines.
3. **The kernel is atomic-bound, and locality is worth 2.9× for free.** Sweeping
   only the hash table's size at $D = 32768$ — same work, same code, different
   address footprint:

   | table | 8 MB | 32 MB (= L2) | 128 MB | 512 MB |
   | --- | --- | --- | --- | --- |
   | in-kernel | **2.74e9** | 2.19e9 | 1.05e9 | 0.962e9 pair/s |

   The 512 MB row is P4 again; the 8 MB row is what the same kernel does when
   the probes stop leaving L2 (32 MB on this device). Nothing else moved the
   number — vectorizing the `sq[]` loads from 8×`LDG.E` to 2×`LDG.E.128` changed
   it by 0.7%, i.e. not at all.

   The catch is that a smaller table means a smaller window cap, so the window
   count and the planner cost rise together: at 8 MB the planner is 91% of the
   wall clock. **So this 2.9× is unclaimable until the planner is fixed, and
   free the moment it is.** That is the argument for doing the planner first
   rather than the more interesting device work.

4. **This particular champion has no $E_4$.** That is a result about one norm and
   nothing more — S-E is a discovery search, so a null answer excludes nothing.
   It does say the machinery runs at the sized problem, which is what G4 was for.

---

## 4. If S-A ever does go to the device

Only justified past $B \approx 10^{13}$. Recorded now because the campaign taught
things that change the kernel design, and they will be forgotten otherwise.

**The shared-memory budget was overestimated by 10×.** The design doc budgets the
ancestor path at $\text{depth} \times D_{\max} \times 8$ B = 20 KB. That is not
how the spine grows: $D$ roughly *doubles* per level, so the packed path costs
$\sum_k D_k \approx 2 D_{\max}$, not $\text{depth} \times D_{\max}$. The worst
real case at the design bound ($5^2\cdot13\cdot17\cdot29\cdot37\cdot41\cdot53\cdot61$,
depth 9) has level sizes $1,2,3,6,12,24,48,96,192$ — **384 points, 3 KB with
`int32` coordinates.** Long chains are the other extreme and no worse: $5^{17}$ is
depth 17 with $D = 9$ and sums to ~81 points. So **3–4 KB per block, not 20**,
which moves the occupancy limit from 2–4 blocks/SM to 12–16 and makes subtree
granularity a much less delicate knob than §4.3 of the design doc assumes.

**One thread per norm, not one thread per point.** `implementation_design.md`
§10.1 measured why: the parent point set averages **1.58 points** across the real
tree (58.5% have one, 31.8% have two). Parallelism across the points of a norm
does not exist. The width is across *norms*, and SIMT supplies it — this is
precisely the axis CPU SIMD cannot reach, because the spine makes consecutive
norms dependent. It is also the strongest remaining argument for a device port at
all.

**The existing meet kernel is obsolete in shape.** `E4MeetDirectKernel` is a
runtime-bounded $O(D^4)$ quadruple loop recomputing `PairKey` inside the inner
loop. The CPU meet since superseded it on three counts, all measured: hoist
$sq[\,]$ once per norm (1.3×), dispatch on $D$ with `constexpr` partition tables
(1.6–2.1×, and ~20× projected at the design bound), and **`meet_direct<5>` and
`<7>` must exist** — $|U| = \lceil\frac12\prod(b_i+1)\rceil$ is odd whenever every
exponent is even, first at $4225 = 5^2\cdot13^2$, and a dispatch with cases 4/6/8
sends every odd-$D$ norm to the slow path. A device meet should be a
transliteration of `src/cpu/meet/meet.h`, which was written against a `HitSink`
concept for exactly this reason.

**Segments map to persistent blocks, not to a static grid.** The campaign's unit
is "all norms whose largest prime factor is $p$" — 4.3 M of them at $10^{10}$,
largest 2.2% of total cost, which is what makes parallelism possible at all
(89.4% for the smallest-prime alternative). But cost falls steeply with segment
index — the first 4,096 segments are 39.4% of the work — so a static
block-per-segment grid is catastrophically imbalanced. Persistent blocks pulling
from a device-side atomic counter is the direct analogue of the CPU's thread pool,
and the adaptive batching (§9.1) is the analogue of the checkpoint interval.

**Keep the false-sharing lesson.** Per-thread campaign state needed
`alignas(64)`: ~470-byte objects adjacent in a vector put one thread's per-norm
counters on the same line as the next thread's spine header. ThreadSanitizer
cannot see this — it is not a race — and it will not see the device analogue
either. On the GPU the same class of bug is bank conflicts in the shared-memory
spine and contention on the append counter; neither shows up in
`compute-sanitizer`, and both need `ncu` sector/conflict counters.

**Precision unchanged.** `uint64` filter on device, `__int128` confirm on host,
GMP certificate on host. Certification is ~25 µs/solution against 8 solutions at
$10^{10}$ — it never belongs on the device.

---

## 5. Existing GPU code: keep, fix, retire

| Component | State | Action |
| --- | --- | --- |
| `gpu/core/` `DeviceBuffer`, `HostBuffer`, `check.h`, `launch.cuh` | sound | **Keep** |
| `gpu/core/device_append_buffer.h` | sound; the "every append bumps the counter" contract is right, and the `__CUDACC__` ODR note is accurate | **Keep** |
| `gpu/meet/` `E4MeetDirectKernel` | superseded on every axis (§4) | **Retire** when a real meet lands; keep the launch/append plumbing |
| `gpu/bigint/` | schoolbook to 1024 bits, one thread per operand pair | **Promote** — S-E needs a fixed 4-limb type built on this convention |
| `gpu/matmul/` | pure capability demo | Leave as sandbox |
| **Pinned-memory pooling** | **[measured]** 1.45 GB/s allocating per transfer against **13.4 GB/s pooled — a 9.2x penalty** (P7). Pageable is 13.1 GB/s, so the win is *not allocating*, not pinning | ~~Unfixed defect~~ **Fixed** — `gpu/core/pinned_pool.h`, one allocation reused across every transfer; its test asserts on `allocations()`, not on a timing. `gpu/champion` uses it and the window planner writes straight into it. `gpu/meet` still allocates per call and is slated for retirement anyway |
| `__launch_bounds__` tuning in `meet_impl.h` | measured, correct, and about to be irrelevant | Re-measure per kernel; the recorded table is for a host-fed kernel |

Two invariants from the design doc that must survive any port: `atomicAdd` takes
exactly `unsigned long long` (on LP64 `std::uint64_t` is `unsigned long`, a
distinct type), and `double` runs at 1/64 rate on sm_89 so every kernel stays
integer-only.

---

## 6. Correctness protocol

Unchanged in substance from `gpu_algorithms.md` §7, but the oracle is now real
rather than aspirational, so the protocol gets sharper:

- **Differential against the CPU catalog.** The $10^{10}$ run is byte-identical
  across thread counts and across eight `SIGKILL`/resume cycles. A device S-A
  must reproduce it exactly after sorting (atomic append makes hit order
  nondeterministic — a test comparing raw output will flake).
- **The audit module applies unchanged.** Known solutions from Bremner's two
  families and `code/old` are checked by enumerating one norm each, so it is
  $O(\text{known})$ work at any bound and any backend. The two above $10^{10}$
  are 33,065,996,717 and 43,937,864,165; F2 $k{=}8 = 2{,}703{,}273{,}548{,}105$
  is the negative control above $10^{12}$.
- **Non-squarefree coverage is the test that finds real holes.** Two of the eight
  catalogued solutions are non-squarefree ($37^2\cdot101\cdot401$ and
  $53\cdot73^2\cdot857$). Branching over each *copy* of a repeated prime is what
  produces them, and a device generator that iterates distinct primes realises
  only exponent 0 or $b$. Norm counts will not catch this; the family self-test
  will.
- **Canonization before the meet, always.** $\nu(\bar z) = \nu(z)$, so on raw
  products $\nu$ is not injective, and $\nu$-injectivity is the sole licence for
  omitting the disjointness test from every kernel. A device pipeline that
  reorders canonization for coalescing silently reintroduces the need for it.
- **For S-D there is no oracle**, which is why §2.1's brute-force check against
  polynomial factorization over $\mathbb{F}_p$ matters and must ship as a test.
- **For S-E there is no oracle either**; the check is a planted collision —
  inject a known equal-key quadruple into a champion norm's point set and require
  recovery, plus a small-$D$ run against `meet_oracle`.

---

## 7. Milestones

Ordered by *information per unit of work*, which is not the order either earlier
draft used.

| | What | Why here | Gate |
| --- | --- | --- | --- |
| **G0** ✅ | ~~Run the profiling queue~~ **Done.** `code/native/bench/` is a built suite with a shared harness and `perf_event_open` counters; results in `drafts/performance_ceiling.md`. It refuted one hypothesis (sorting networks in `Spine::push`) and promoted another (`meet_hashed` is 82% sort) | Nothing in this repo had a real profile. §10.1 found the cost sitting in a three-element `std::sort` — intuition here is demonstrably unreliable | none |
| **G1** ✅ | ~~Settle the S-D local model~~ **Done.** The ring is $\mathbb{Z}/p^k$ and the missing parameter is a cap on node valuations; the two models §2.2 weighs are its two extreme cases at $k=1$. `drafts/local_model_plan.md`, `paper.md` §6.6 + Prop. 4.3, `code/native/src/cpu/tower/` | Blocked $\sigma_p$ and the wheels. It also produced arithmetic: computed forced divisors that exceed [CMSV]'s $C'_{16}$ at $p = 2, 3, 11$, and are sharp at four of six small primes at depth 4 | G0 not required |
| **G2** ◐ | S-D DP on CPU **with the brute-force fixture** — done; **to the precision floor** — not | The fixture ships (exhaustive tree-walk oracle, exact agreement at every cap). What is missing is reach: the floor is $k \ge 15$–$20$ at $p=3$, the current $O(p^{3k}/64)$ costs 63 s at $3^8$, and $3^9$ takes hours. Needs the §5.1 rewrite | G1 ✅ |
| **G3** ✅ | ~~Fixed 4-limb integer type~~ **Done.** `code/native/src/core/uint.h` — `UInt<N>` over **32-bit** limbs (the convention `cpu/bigint` and `gpu/bigint` already state), with `add.cc.u32`/`addc.cc.u32` carry chains on the device; `gpu/core/pinned_pool.h` for the staging | Both prerequisites for S-E. The limb width was not a free choice: **P8** measures 7 SASS instructions per 256-bit add for the 32-bit PTX form against ~25 for a 4×`uint64` C one, because with 64-bit limbs the compiler stops using the carry flag and rebuilds it with ISETP+SEL. 2.2× at equal clock, 0.54 s across all of S-E | none |
| **G4** ✅ | ~~**S-E**: value-windowed two-pointer meet~~ **Done.** `cpu/champion/` (reference, exact 256-bit keys) and `gpu/champion/` (the hash-insert kernel), cross-checked in `apps/champion_demo.cc` | The clearest GPU win, and it held: §3.1's windowing is implemented as an adaptive bisection with the exactly-once property tested against brute force, and the sized champion runs on the device. See §3.4 for what it measured | G3 |
| **G5** | S-D on device, only if $\sigma_p$ converges slowly enough to need $p > 10^3$ | Conditional by construction | G2 |
| **G6** | S-C bottom-lift, sized on CPU before any kernel | Currently the only search with *no* cost estimate at all | catalog |
| **G7** | Device S-A — K2/K3 first, against the existing catalog | Only if $B$ moves past $10^{13}$. P1 puts the device meet's ceiling at 7.34e11 quad-partitions/s, ~46× sixteen CPU threads — but the design bound's whole bucketed meet is only 1.48e12 partitions, i.e. **~2 s of device arithmetic**. Generation is the cost, not the meet | G0 ✅ |

G0–G2 involve no CUDA at all, and that is the finding, not an oversight.
**G0, G1, G3 and G4 are now complete**; G2 is half done and blocked on an
algorithmic rewrite rather than on hardware. G5 and G7 remain conditional by
construction, and G6 still has no cost estimate at all.

The next unblocked item is **G2's second half**, and it is worth naming what it
is *not*: not a CUDA task, not a threading task, but the orbit-representative
rewrite of `local_model_plan.md` §5.1. The ordering matters because the factors
multiply in the wrong order otherwise — the rewrite is $m/64$ (≈300× at $3^9$),
threads are ≤16×, and the device is maybe another 10×, but the rewrite also
*deletes* the irregular structure the other two would have had to parallelise.

Four items the S-E build added to the queue. The first three are sized rather
than guessed and are in this order because each one gates the next; the fourth
is not engineering at all:

- **G8 — get the planner off the critical path.** Measured at exactly 50% of the
  sized run (§3.4), with the device idle throughout. A second stream planning
  window $N{+}1$ while $N$ runs is ~2× for a few dozen lines; a device-side
  planner is more, and is what makes small windows cheap.
- **G9 — then shrink the table into L2.** A measured 2.9× that costs no new
  code, and is unaffordable today only because it multiplies the window count
  (§3.4). G8 converts it from a losing trade into a free one.
- **G10 — then, if it still pays, partition the hash.** Shared-memory
  write-combining buckets flushed to global when full, re-bucketed in a second
  kernel where the fan-out exceeds what shared memory holds — the standard
  multi-pass radix-partitioned join. It attacks both halves at once: coalesced
  flushes and one counter atomic per line instead of per key, *and* a window
  bounded by a key buffer rather than by table capacity, which cuts the window
  count ~5× and so the planner with it. Napkin: ~4 GB of streaming per 1.3e8
  keys → ~4e9 key/s, against 2.7e9 for G9's L2-resident table. **The honest
  comparison is against G9, not against today**, which makes it perhaps 1.5×
  for real work — worth doing after G8/G9 and not before, and P3's CUB radix
  sort at 8.4e8 key/s is a caution that multi-pass streaming on this device is
  not automatically ahead of direct atomics.
- **Which champion?** Open question 5 is now the operative one for S-E rather
  than a hypothetical: the machinery runs, so the next result depends entirely
  on which norm it is pointed at. `smallest_split_primes` answers only the
  obvious version (maximise $D$ per unit $L_1$).

---

## 8. Open questions

1. ~~**What is the right local model for $\sigma_p$?**~~ **Answered (§2.2), and
   the answer was that the question was posed over the wrong ring.** Both models
   are the $k=1$ truncation of a $p$-adic quantity: a solution is a
   $\mathbb{Z}_p$-point of the *nondegenerate* locus, nondegenerate means every
   node radical has **finite** valuation, and mod $p$ one can see zero and
   nonzero but not finite. Over $\mathbb{Z}/p^k$ with a cap $c$ on node
   valuations, $c=0$ *is* the distinct-root model and $c=k$ *is* the
   multiplicity model. What replaces the question, in order:
   (a) the nonexistence branch needs $V_p(k,d) \to \infty$, a limit no finite
   $k$ settles — the missing ingredient is eventual periodicity or Poincaré
   rationality, and it is still mathematics;
   (b) $\sigma_p$ additionally needs a coordinate correction ($L$-count vs. root
   count) that is not yet done;
   (c) the wheels are available now and are proven congruences.
   **The new highest-value open question is (a).**
2. ~~**Does the value-windowed two-pointer (§3.1) extend to the exhaustive
   regime?**~~ **[measured] Answered, and the guess in this slot was wrong.** The
   sort is *not* already cheap: `meet_hashed` is **81.9% sort**, carrying 86
   branch misses per norm and collapsing IPC from 5.55 to 1.53
   (`performance_ceiling.md` §2.2). Since $D>8$ is 3% of norms but 74% of
   direct-meet work, and the meet is 57% of the inner loop at $10^9$ and rising,
   **a sort-free `meet_hashed` is now the highest-value open optimization in the
   project** — the same structural move as §3.1, one level down.
3. **Is A4′'s index-calculus framing competitive?** Unchanged from
   `gpu_algorithms.md` §6: collisions are equalities of factorizations over a
   prime set confined to $\equiv\pm1 \bmod 8$ by Theorem 5.3. Still the only route
   past the champion regime's quadratic wall, still open.
4. **What does S-C actually cost?** Divisor enumeration of $2\mu r'_k$ per root,
   joined on shared $L$, with the parity wheel killing $\mu = 1,2,3$. No estimate
   exists in any draft. It could be the cheapest search or the most expensive.
5. **Which champion norms?** Maximise $D$ per unit work — primorial-style products
   of the smallest split primes, or tuned to the §4 ladder congruences? Unchanged,
   and now sized: 20 distinct primes gives $D = 2^{19}$ at $L_1 \approx
   2.5\times10^{36}$.
