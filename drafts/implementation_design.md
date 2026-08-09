# Implementation Design — CPU and GPU Enumeration

**Scope.** How the searches of `paper.md` §6 and the algorithms of
`gpu_algorithms.md` §4 become a program. Part I fixes *what* we are building and
at which working point; Part II says *where the time goes*, quantitatively, so
the structure in Part III is forced by measurement rather than taste; Part IV is
the list of things that will bite.

Every number below is either measured on this machine or an exact sieve count.
Extrapolations are labelled as such and carry their growth factor.

Hardware of record: RTX 4060 Laptop (sm_89, 8 GB, 3105 MHz, **PCIe 4.0 ×8**),
CUDA 12.9, g++-13. Single-core CPU rates are from `RelWithDebInfo`.

---

# Part I — What we are implementing

## 1. The algorithm inventory

Five searches, in dependency order. S-A is the spine: everything else consumes
its catalog or prunes its input.

### S-A — Per-norm $E_4$ enumeration (the exhaustive campaign)

For each admissible $L_1 \le B$ build $U(L_1)$, take $\nu(z) = \mathrm{Re}(z^2)$,
find 4-subsets splitting as $\{z_1,z_2\}\,|\,\{z_3,z_4\}$ with
$\nu_1^2+\nu_2^2 = \nu_3^2+\nu_4^2$.

| Working point | Value |
| --- | --- |
| Target bound $B$ | $10^{12}$ (design), $10^{10}$ (first campaign) |
| Factor-base cutoff | $P_{\max} = B/65$ — exactly tight, see §4.1 |
| Admissible norms at $B{=}10^{12}$ | $2.25\times10^{10}$ (extrap., $11.1\times$/decade) |
| Points | $1.26\times10^{11}$ (extrap., $11.5\times$/decade) |
| Meet work, direct | $2.97\times10^{13}$ partitions (extrap., **$20.5\times$/decade**) |
| Meet work, bucketed | $1.48\times10^{12}$ (extrap., $13.3\times$/decade) |
| $D_{\max}$ | 192, at $5^2\cdot13\cdot17\cdot29\cdot37\cdot41\cdot53\cdot61$ |
| DFS depth | $\le 11$ measured at $10^8$; $\le 13$ expected at $10^{12}$ |
| Arithmetic | `uint64` filter + `__int128` confirm (§3 of `gpu_algorithms.md`) |

### S-B — Top-collision join

Hash-join the S-A catalog on $(L_1,L_2,L_3)$; a collision with distinct $L_4$ is
an $E_5$. Working point: catalog size, not arithmetic — radix-partition by a
hash of the key so each partition fits in RAM. Depends on S-A completing.

### S-C — Bottom-lift join

Per catalogued $E_4$ with roots $r'_0..r'_7$, seek $(\mu, L)$ with $L \pm \mu r'_k$
both squares for all $k$. Working point: divisor enumeration of $2\mu r'_k$ per
$k$, joined on shared $L$. Parity wheel kills $\mu = 1,2,3$ outright (odd roots
force $\mu \equiv 0 \bmod 4$).

### S-D — $\mathbb{F}_p$ tower sieve

DP over tower levels, $O(p^2)$ per level — **not** $p^5$ brute force. Working
point: $p$ up to $\sim10^5$–$10^6$. Outputs local densities $\sigma_p$ for the
§8 heuristic, sieve wheels for S-A, and the jackpot branch (one prime with zero
depth-5 towers is an unconditional nonexistence theorem). **Independent of S-A**
— buildable and runnable in parallel.

### S-E — Champion-regime meet

Single norms with $D \sim 2^{19}$, where $D^2/2 \approx 1.4\times10^{11}$ keys
cannot be materialized. Windowed hashing: keep only keys whose hash lands in a
residue window of width $1/R$, repeat $R$ times. At $R = 10^3$ each pass retains
$\sim1.4\times10^8$ keys ($\sim$1 GB), total $\sim1.4\times10^{14}$ key
evaluations.

### Support: A5 confirmation and tower recovery

On every filter hit: recompute exactly, peel $(L_1,L_2,L_3,\sqrt{L_4})$, assert
the degree-16 identity, the ladder divisibilities, and primitivity. Rare, so its
cost is irrelevant and its correctness is everything. Measured at **~25 µs per
solution** in GMP, against ~10⁻⁷ s per norm for the search — so certification is
free at any bound where solutions stay as sparse as they are (2 below $10^6$).

The peel is exact and unconditional, which is worth stating because it means the
certificate cannot *discover* that a hit was not a solution: $\nu^2 + u^2 =
L_1^2$ makes the meet's $\nu$ condition and the tower's $u$ condition the same
equation, so $L_2$ exists for every hit and $L_3$, $\sqrt{L_4}$ are then
determined. What it can discover is that *this code* computed the wrong tower —
which is exactly what a paper artifact needs it to catch.

## 2. Workload shape (exact sieve counts)

| $B$ | factor base | norms | points | pairs | quad-partitions | $D_{\max}$ | $D{=}4$ | $D\le8$ | $D$ odd | $D>8$ work |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| $10^6$ | 891 | 11,083 | 50,617 | 99,354 | 347,919 | 16 | 79.4% | 98.84% | 1.89% | 35.3% |
| $10^7$ | 7,072 | 131,437 | 623,030 | 1,339,509 | 8,295,267 | 24 | 77.3% | 97.88% | 1.24% | 59.5% |
| $10^8$ | 58,381 | 1,464,114 | 7,177,311 | 16,790,533 | 169,730,031 | 36 | 75.5% | 97.03% | 0.87% | 74.2% |

"Quad-partitions" is $\sum_{\text{norms}} \binom{D}{4}\cdot 3$ — the actual work
of the direct meet, and the metric that matters.

> **Correction.** An earlier version of this table gave $D\le8$ as
> 97.4 / 97.0 / 96.5%. Those are the shares of $D \in \{4,6,8\}$: the script
> behind them assumed $D$ is always even. It is not.
> $|U| = \lceil\tfrac12\prod(b_i+1)\rceil$, so any norm with all-even exponents
> has an odd point count — $4225 = 5^2\cdot13^2$ is the smallest, at $D = 5$.
> The `$D$ odd` column above is what was being silently dropped. The consequence
> for §7.1 is concrete: a dispatch with cases 4, 6, 8 sends every odd-$D$ norm
> to the large-$D$ path, so `meet_direct<5>` and `<7>` exist in the shipped
> code.

### 2.1 Where the points come from, and why canonization is not optional

The multiplicative construction of $U(L_1)$ generates, for a factorization with
$k$ prime *copies*, exactly $2^k$ products — one per choice of
$\mathfrak{q}_i$ or $\bar{\mathfrak{q}}_i$ at each copy. Those are not $2^k$
distinct points. Canonizing to the first octant ($m \ge n \ge 0$) collapses each
orbit under the units and conjugation, leaving

$$|U(L_1)| = \left\lceil \tfrac{1}{2}\textstyle\prod_i (b_i+1) \right\rceil .$$

The halving is conjugation: the construction already emits one representative
per unit class, and $z$ and $\bar z$ appear as the complementary exponent
choices $c \leftrightarrow b-c$. The $\lceil\cdot\rceil$ covers the
self-conjugate point, which exists only when every $b_i$ is even. Verified:

| $L_1$ | factorization | raw products $2^k$ | $|U|$ |
| --- | --- | --- | --- |
| 1105 | $5\cdot13\cdot17$ | 8 | **4** |
| 32045 | $5\cdot13\cdot17\cdot29$ | 16 | 8 |
| 5525 | $5^2\cdot13\cdot17$ | 16 | 6 |

So **three distinct primes yield four independent points** — the minimum for an
$E_4$, which is exactly why $5\cdot13\cdot17 = 1105$ is the smallest admissible
norm and why $D = 4$ is the modal case.

**Now the part that is easy to get wrong.** It is tempting to think canonization
is mere deduplication that $\nu$ makes unnecessary — after all the key uses
$\nu^2$, so it is blind to the sign flip $\nu(iz) = -\nu(z)$. That reasoning is
wrong, and the failure is not a performance loss but a correctness one:

$$\nu(\bar z) = \nu(z).$$

$\nu$ is *invariant* under conjugation. On the raw products, therefore, **$\nu$
is not injective** — $z$ and $\bar z$ carry the same value. Measured on all
three cases above: injective on the canonized set, never on the raw one.

And $\nu$-injectivity is precisely what licenses omitting the disjointness test
from every meet, on CPU and GPU alike. Without canonization, two pairs with
equal keys could share a point: $\{z_1, z_2\}$ and $\{z_1, \bar z_2\}$ have equal
keys and overlap in $z_1$. The kernel would report it as a solution, and it is
not one.

> **Canonization is what makes $\nu$ injective, and $\nu$-injectivity is what
> removes the disjointness test.** Anything that reorders the pipeline so points
> reach the meet uncanonized silently reintroduces the need for that test.

It is also a $2^k/|U| \approx 2\times$ reduction in $D$, hence a **~16×**
reduction in $O(D^4)$ meet work — but that is the lesser reason.

**Two facts set the entire design.**

**(a) $D$ stays tiny.** Three quarters of admissible norms have exactly four
points. $D \le 8$ covers ~96.5%. At $D=4$ there are exactly three disjoint
pair-partitions — three comparisons, no hash table, no shared memory, no
disjointness test. This is a throughput problem over $\sim10^{10}$ near-identical
micro-tasks, not a compute-bound one.

**(b) The 3% tail dominates the cost, and increasingly so.**

| $B$ | $D>8$: share of norms | $D>8$: share of direct-meet work |
| --- | --- | --- |
| $10^6$ | 1.2% | 35.3% |
| $10^7$ | 2.1% | 59.5% |
| $10^8$ | 3.0% | **74.2%** |

$\binom{D}{4}$ growth outruns the thinning of the tail. This is why direct-meet
work grows at $20.5\times$/decade while norms grow at $11.1\times$ — and it is
the quantitative case for $D$-bucketing:

| | $10^6$ | $10^7$ | $10^8$ | $10^{12}$ (extrap.) |
| --- | --- | --- | --- | --- |
| direct | 347,919 | 8,295,267 | $1.70\times10^8$ | $2.97\times10^{13}$ |
| hash for $D>8$ | 231,825 | 3,540,276 | $4.71\times10^7$ | $1.48\times10^{12}$ |
| **speedup** | 1.5× | 2.3× | 3.6× | **~20×** |

Bucketing is not a micro-optimization deferred to later. It is the difference
between an algorithm that grows like $B^{1.3}$ and one that grows like $B^{1.1}$.

---

# Part II — Where the time goes

## 3. CPU critical path

Measured single-core, at $B = 2\times10^6$ (23,576 norms, 108,909 points,
925,452 quad-partitions):

| Stage | Time | Unit rate |
| --- | --- | --- |
| Factor base: sieve + Miller–Rabin + Cornacchia | 0.85 ms | 0.5 µs / prime (1,649 primes) |
| Point-set construction | 3.63 ms | **154 ns / norm**, 33 ns / point |
| Reference meet | 1.85 ms | **2.00 ns / quad-partition** |

Extrapolating to the design bound at these rates:

| $B$ | factor base | point sets | meet (direct) | meet (bucketed) |
| --- | --- | --- | --- | --- |
| $10^{10}$ | ~14 s | 0.05 core-h | 0.04 core-h | **0.005 core-h** |
| $10^{12}$ | ~164 s | 0.96 core-h | 16.5 core-h | **0.95 core-h** |

**The headline, and it should change the plan:** with bucketing, the exhaustive
$E_4$ campaign to $B = 10^{12}$ is on the order of **two core-hours of
arithmetic** — call it 10 minutes on 16 threads. It is not a supercomputing
problem.

One caveat, stated because the estimate is load-bearing: the 154 ns/norm rate
will degrade. It is dominated by two `std::vector` allocations per `point_set`
call, and at $10^{12}$ average $D$ is ~5.3 rather than 4.6, with far worse cache
locality across a $10^{11}$-point stream. Assume 2–3× pessimism.

## 3.1 One-time work vs. per-campaign vs. per-norm

The distinction matters more than any single rate, because it decides what is
allowed to be expensive. Three tiers:

| Tier | Runs | Budget |
| --- | --- | --- |
| **One-time (OTW)** | once, ever, per bound; result is an on-disk artifact | minutes are free |
| **Per-campaign** | once per run | seconds |
| **Per-norm** | $2\times10^{10}$ times | nanoseconds |

**The factor base is OTW, and that changes its status entirely.** An earlier
draft of this document called it "the real wall". That was wrong, on both
counts, and the measurements are worth recording because they invert the
conclusion.

Measured, with a segmented sieve plus Cornacchia, single-threaded:

| $B$ | $P_{\max}$ | split primes | sieve | Cornacchia | total | **resident** | artifact |
| --- | --- | --- | --- | --- | --- | --- | --- |
| $10^8$ | $1.54\times10^6$ | 58,433 | 0.01 s | 0.04 s | **0.05 s** | 0.1 MB | 0.4 MB |
| $10^9$ | $1.54\times10^7$ | 497,241 | 0.07 s | 0.39 s | **0.46 s** | 0.1 MB | 3.8 MB |
| $10^{10}$ | $1.54\times10^8$ | 4,327,911 | 0.73 s | 3.86 s | **4.59 s** | 0.1 MB | 33 MB |
| $10^{11}$ | $1.54\times10^9$ | 38,308,127 | 7.33 s | 38.74 s | **46.1 s** | 0.2 MB | 292 MB |
| $10^{12}$ | $1.54\times10^{10}$ | 343,618,258 | 74.8 s | 462.1 s | **536.9 s** | 0.2 MB | 2.6 GB |

The last row is **built and verified**, not extrapolated — `data/factorbase_1e12.bin`,
343,618,258 records, 8.9 minutes single-threaded, 0.2 MB resident throughout.

Two things kill the "wall":

1. **Segmentation makes resident memory constant.** 0.1–0.2 MB regardless of
   $P_{\max}$ — the segment plus the base primes to $\sqrt{P_{\max}}$ (1,481
   primes at $B{=}10^{10}$, 4,134 at $10^{11}$). The 1.9 GB figure was an
   artifact of the *current monolithic* `std::vector<bool>`, which is a defect
   in `split_primes_up_to` (**C8**), not a property of the problem.

2. **The right encoding is $(a,b)$, not $(p,a,b)$.** $p = a^2 + b^2$ is
   recoverable, so the prime never needs storing. At $B = 10^{12}$,
   $a, b < 1.24\times10^5$ — 17 bits each — so a pair of `uint32` is 8 B/prime
   and the whole base is **2.7 GB**, not 5.2. On disk, `mmap`'d, streamed
   sequentially: it never needs to be resident, and the page cache handles the
   access pattern for free.

So: build it once, ship it as an artifact, `mmap` it. ~8 minutes single-core, and
the work is embarrassingly parallel over segments — under a minute on 16 cores.

The residual truths, both narrower than "wall":

- `split_primes_up_to` as written really will not run at $B \gtrsim 10^{10}$.
  That is a code defect to fix, worth ~1 day, not an algorithmic barrier.
- 2.7 GB is more than an 8 GB GPU should dedicate, so the device streams
  segments of the artifact rather than holding it. At PCIe rates that is 0.2 s
  for the whole base — irrelevant.

**Consequence for the GPU decision.** The exhaustive regime does not need the
GPU. Its value is S-D (the $\mathbb{F}_p$ sieve, enormous intensity, independent
of S-A) and S-E (champion norms). This is `gpu_algorithms.md` §5 option A, now
with numbers behind it.

## 4. GPU critical path

### 4.1 What crosses PCIe: only the factor base

The device generates its own norms and point sets. Nothing else makes sense
once the numbers are in front of you, so state the traffic budget first and let
the pipeline follow from it.

| Direction | Payload | Volume at $B = 10^{12}$ |
| --- | --- | --- |
| H→D, once per segment | decomposed primes $(a,b)$, `uint32` pairs | **2.7 GB total**, streamed |
| H→D, once | cofactor `SmallTable` (DFS tree, ~88k nodes) | ~2 MB |
| D→H | confirmed candidates | ~0 |

That is it. The factor-base artifact of §3.1 *is* the input: 8 bytes per prime,
sequential, `mmap`'d on the host and copied segment by segment. At PCIe rates
the entire base crosses in **0.2 s** — for a campaign whose arithmetic is hours.
Traffic stops being a design constraint at all.

Contrast the naive hybrid, which is what the current demo does and why it exists
as a demonstration rather than a plan. Shipping $\nu$ values instead of primes,
at the modal $D = 4$ norm:

```
payload   4 nu values x 8 B                    =  32 B
work      3 partitions x (2 mul + 1 add) + cmp ~  15 int64 ops

intensity 15 / 32                              = 0.47 ops/byte
transfer  32 B / 12 GB/s                       = 2.667 ns
compute   15 / 1e12 ops/s                      = 0.015 ns
                                                 ---------
                                                 178x transfer-bound
```

Aggregate: $1.26\times10^{11}$ points × 8 B = **1008 GB**, ~84 s of pure PCIe,
to accelerate ~1 core-hour of work. The measured demo agrees at 0.25 ops/byte.
**373× more traffic than sending primes.** The naive hybrid is not a slower
design; it is a wrong one, and no amount of tuning reaches the right one from
there.

### 4.2 The kernel pipeline, and the milestone order

Multi-kernel, with each stage writing to global memory and the next reading it.
Fusion is an *optimization applied later*, not the starting point — the staged
version is debuggable, each kernel is independently testable against the CPU,
and the intermediate arrays are exactly what a cross-check needs to look at.

```
  [host]  mmap factor-base artifact, pick a segment of primes
     |    H->D: (a,b) pairs for the segment              ~8 B/prime
     v
  K1  expand cofactor subtrees -> norm descriptors
      one block per subtree; emits (prime indices, exponents, D)
      output: norm list, bucketed by D (free: D(cq) = 2D(c))
     |
     v
  K2  build point sets
      one warp (or block) per norm; U(cq) = U(c)q  u  U(c)q-bar
      output: nu values in CSR layout, per D-bucket
     |
     v
  K3  the meet, one kernel instantiation per D-bucket
      meet_direct<4>, meet_direct<6>, meet_direct<8>, meet_hashed
      output: candidates, atomically appended
     |
     v
  [host]  confirm in __int128, primitivity, tower recovery
```

**Milestone order**, and the first two are worth shipping on their own:

| Milestone | What it is | Why this order |
| --- | --- | --- |
| **M1** | K2 + K3 as two kernels, norms still from the CPU | Isolates the point-set kernel — the new, untested part — while the meet is already validated. `nu` lands in global memory where the existing CPU oracle can read it verbatim. |
| **M2** | K1 on device; nothing but primes crosses | Removes the transfer entirely. Now the whole pipeline is device-resident and the §4.1 budget is real. |
| **M3** | Fuse K2+K3 for the $D \le 8$ buckets | Point sets for 96.5% of norms fit in registers, so materializing `nu` to global is pure waste — but only *after* M1 has proved the two halves agree separately. |
| **M4** | Fuse K1 into the block: subtree walk in shared memory | The full §4.3 design. Highest payoff, hardest to debug, so last. |

The bucketing is what makes K3 clean: because $D(cq) = 2D(c)$, the bucket is
known before a point set is built, so K1 can emit *already partitioned* norm
lists and K3 launches one homogeneous grid per bucket. No stream compaction, no
divergence, and each grid is a different template instantiation.

### 4.3 Where each stage lives, and why

The unusual thing about this workload is that the memory hierarchy is used
*backwards* relative to normal GPU practice. Inputs are tiny and streamed;
outputs are essentially nil; the pressure is entirely on **intermediates**. The
usual instinct — optimize output bandwidth, coalesce the writes — optimizes the
one thing that does not matter here.

| Stage | Lives in | Size | Reasoning |
| --- | --- | --- | --- |
| Factor-base segment $(a,b)$ | global, read-once | 8 B/prime, ~10–100 MB/segment | Streamed sequentially and never revisited. Coalesced by construction; no reuse to exploit, so L2 is irrelevant. |
| `SmallTable` cofactor tree | global, read-mostly | ~2 MB | Every block reads it, none writes it. Small enough to stay hot in L2 across the whole kernel — this is the one input with real reuse. |
| Subtree walk state | **shared**, per block | ≤ 20 KB | A block owns one subtree and walks it emitting thousands of norms. The ancestor path is depth ≤ 13 × $D$ ≤ 192 × 8 B = 20 KB, against 48–96 KB/SM — so 2–4 blocks/SM by shared memory alone, and far less for typical subtrees. |
| Point set $U$, built incrementally | **shared** | $D$ × 8 B | A child is its parent doubled by one prime: $U(cq) = U(c)\mathfrak{q} \cup U(c)\bar{\mathfrak{q}}$, one Gaussian multiply per point. Rebuilding from scratch instead costs a factor of the tree's branching. |
| $\nu$, squares, keys | **registers** | $D$ × 8 B | 32 B at $D=4$. This is where 96.5% of all arithmetic happens and it never touches memory at all. |
| Candidates | global, atomic append | ~0 | Hits are so rare that the atomic is uncontended and the write bandwidth is unmeasurable. Do not design for it. |

The two consequences worth internalizing:

1. **Occupancy is bounded by shared memory, not registers**, once the subtree
   walk moves on-device. The 20 KB path is the binding constraint, so subtree
   granularity — how deep a subtree one block owns — is the real tuning knob,
   not threads-per-block.
2. **The arithmetic is free.** At $D = 4$ the meet is 4 squares and 3 compares
   against a point-set build of 3 Gaussian multiplies. Every kernel here is
   latency- and scheduling-bound, never throughput-bound, which is why `double`
   at 1/64 rate is irrelevant and why integer instruction counts barely move the
   needle.

### 4.4 Occupancy tuning — deferred

`E4MeetDirectKernel` currently sits at 56 registers, 0 spills, 9 blocks/SM, 75%
occupancy on sm_89, with `__launch_bounds__(128, 9)`. Recorded in
`meet_impl.h` with the measured table behind it.

This is **not** design input, and the numbers will change the moment the kernel
stops being fed from the host. Treat it as a placeholder for a later autotuning
pass over (threads/block, subtree depth, bucket) once M2 lands. The one durable
lesson from measuring it: `__launch_bounds__(N, 1)` is *not* "unconstrained" —
it promises one resident block, so nvcc spends more registers and occupancy
drops. Always read `-Xptxas -v` rather than reasoning about it.

---

# Part III — Software structure

The structure follows from Part II:

- $D$ is a **compile-time** quantity for 96.5% of work → the meet must be a
  template on $D$ with fully unrolled partition tables, not a runtime loop.
- Generation and meet **start as separate stages and fuse only at M3/M4**
  (§4.2), so the algorithm must be written once against an abstract sink — a
  sink that writes `nu` to global memory in the staged version and one that
  keeps it in registers in the fused version, with the *same* meet source
  compiled against both. Getting this seam right early is what makes fusion a
  later switch rather than a rewrite.
- Precision is a **policy**, not a type choice scattered through the code.
- The DFS spine carries incremental state → an explicit stack object, not
  recursion with rebuild.

## 5. Layering

```
core/         precision policy, strong types, concepts           (header-only)
arith/        GaussT<T>, Cornacchia, Miller-Rabin, isqrt          (backend-agnostic)
factorbase/   segmented sieve + decomposition, as a range
generate/     cofactor DFS, incremental point sets                (the spine)
meet/         MeetDirect<D>, MeetHashed, bucket dispatch
confirm/      exact re-check, tower peeling, primitivity
catalog/      hit sink, serialization, resumable checkpoints
backend/      HostBackend, DeviceBackend                          (the seam)
```

`arith/`, `meet/` and `confirm/` are written once and instantiated for both
backends. `generate/` has two implementations sharing the state machine.

## 6. Core abstractions

### 6.1 Strong types — make the filter/exact distinction unforgeable

The single most dangerous confusion in this program is using a truncated key
where an exact one is required. Types should make it impossible rather than
document it.

Two questions this raises, both worth answering explicitly because the answers
are what make the scheme coherent rather than decorative:

**Why is `Nu` not `TruncNu`?** Because $\nu$ never truncates. $\nu(z) = m^2-n^2$
is bounded by $L_1 \le 10^{12}$ — 40 bits in a 63-bit signed type, exact with
23 bits to spare, at *every* bound this project will run. Only the **key**
$\nu_a^2+\nu_b^2$ overflows, at 81 bits. So the truncation lives in exactly one
place and the name should say so; calling the input truncated too would spread
the warning until it stopped meaning anything.

**Why `__int128` and not `mpz_class` for the exact layer?** Not performance —
confirmation is rare enough that a heap allocation would be affordable. The
reason is that a *fixed* width can be checked at compile time and an arbitrary
one cannot:

```cpp
static_assert(static_cast<exact_key_t>(kMaxNorm) * kMaxNorm * 2 > 0,
              "exact_key_t overflows on 2*kMaxNorm^2");
```

That assertion fires the day someone raises `kMaxNorm` past $9.2\times10^{18}$.
With `mpz_class` there is nothing to assert — the bound becomes untracked, and
"it is arbitrary precision so it must be fine" is a belief rather than a
checked property. Since the whole point of the layer is to be the *arbiter*, a
machine-checked bound is worth more than unbounded range we do not need.

Where the bound genuinely cannot be fixed, GMP is mandatory and there is no
argument: the layer-3 PTE identity needs 311 bits and a size-32 identity 638.
The rule is therefore **fixed width with a static_assert wherever a bound
exists; GMP only where one does not.**

```cpp
template <class Tag, class Rep> struct Quantity {
  Rep v{};
  constexpr explicit Quantity(Rep r = {}) : v(r) {}
  friend constexpr auto operator<=>(Quantity, Quantity) = default;
};

struct NuTag; struct TruncKeyTag; struct ExactKeyTag; struct NormTag;

using Nu       = Quantity<NuTag,       std::int64_t>;
using TruncKey = Quantity<TruncKeyTag, std::uint64_t>;   // mod 2^64, a FILTER
using ExactKey = Quantity<ExactKeyTag, __int128>;        // the ARBITER
using L1       = Quantity<NormTag,     std::int64_t>;
```

`TruncKey` and `ExactKey` share no conversion. The only bridge is an explicit
`confirm(candidate) -> bool` that recomputes in `ExactKey`. A kernel physically
cannot return something a caller mistakes for a result.

### 6.2 Precision as a policy

```cpp
template <class P> concept PrecisionPolicy = requires {
  typename P::coord;  typename P::norm;
  typename P::trunc_key;  typename P::exact_key;
  { P::max_norm } -> std::convertible_to<std::int64_t>;
  { P::exact_key_bound } -> std::convertible_to<std::int64_t>;
};

struct Campaign1e12 {                 // the shipping policy
  using coord = std::int32_t;         // |m| < 1e6: 20 bits
  using norm  = std::int64_t;
  using trunc_key = std::uint64_t;
  using exact_key = __int128;
  static constexpr std::int64_t max_norm = 1'000'000'000'000;
  static constexpr std::int64_t exact_key_bound = 3'037'000'499;
};
```

`coord = int32_t` halves point-set bytes, which matters exactly once
coordinates live in shared memory (§4.3) — the policy is where that flips.

### 6.3 Concepts at the seams

Fair warning on ceremony: §6.1 and §6.2 do overlap, and only one of them earns
its keep universally. **The precision policy is load-bearing** — it is a single
compile-time switch that retargets the whole program from `Campaign1e12` to a
test-scale or a narrow-coordinate variant, and §6.3's backends both key off it.
**The strong types are optional**, and the honest scope for them is narrow: wrap
`TruncKey` and `ExactKey`, because confusing those two is the one mistake that
silently produces wrong answers, and leave `Nu`, `L1` and the rest as plain
integers. Wrapping every quantity buys conversions and noise, not safety.

If it feels like too much machinery when it comes to be written, cut the
`Quantity` template to those two types and keep the policy. The concepts below
are a third mechanism and the lightest of the three — they exist only so the
meet compiles unchanged against a host `std::vector` sink and a device atomic
append.

```cpp
// As shipped: a plain callable, not an object with .emit(). A lambda then works
// with no adapter, which is what the host sink actually is.
template <class S> concept HitSink = std::invocable<S&, const E4Sol&>;

template <class B> concept Backend = requires {
  typename B::template buffer<int>;
  requires HitSink<typename B::sink>;
  { B::name } -> std::convertible_to<std::string_view>;
};
```

The meet is written once against `HitSink`. On the host the sink is a
`std::vector` push; on the device it is `DeviceAppendView::push`. Same source.

## 7. The meet: compile-time partition tables

The core metaprogramming payoff. For $D \le 8$ the set of disjoint pair-pairs is
fixed and small (3 / 45 / 210), so generate it at compile time and emit a fully
unrolled, register-only, branch-free kernel.

```cpp
// Number of ways to split D indices into two disjoint unordered pairs.
constexpr std::size_t partition_count(int D) {
  return static_cast<std::size_t>(D) * (D - 1) * (D - 2) * (D - 3) / 8;
}

struct Partition { std::uint8_t a, b, c, d; };

// {a,b} | {c,d} with a<b, c<d, a<c -- each split generated exactly once.
template <int D>
constexpr std::array<Partition, partition_count(D)> partitions() {
  std::array<Partition, partition_count(D)> out{};
  std::size_t n = 0;
  for (int a = 0; a < D; ++a)
    for (int b = a + 1; b < D; ++b)
      for (int c = a + 1; c < D; ++c) {
        if (c == b) continue;
        for (int d = c + 1; d < D; ++d) {
          if (d == b) continue;
          out[n++] = {std::uint8_t(a), std::uint8_t(b),
                      std::uint8_t(c), std::uint8_t(d)};
        }
      }
  return out;                        // static_assert(n == partition_count(D))
}
```

*This generator is not a sketch — it was compiled and checked against a
brute-force enumeration of pair-partitions for $D = 4,5,6,7,8,10$: exact counts
(3 / 15 / 45 / 105 / 210 / 630), no duplicates, no repeated indices within an
entry, set-equal to truth. The `a` loop taking the global minimum is what makes
each partition appear exactly once.*

```cpp

template <int D, MeetKey Key, HitSink S>
void meet_direct(const norm_t* nu, S&& sink) {
  Key sq[D];
  for (int i = 0; i < D; ++i) {
    const auto v = static_cast<Key>(nu[i]);
    sq[i] = v * v;                                               // A8: hoisted
  }

  constexpr auto table = partitions<D>();
  for (const Partition& t : table) {
    if (sq[t.a] + sq[t.b] == sq[t.c] + sq[t.d])                  // 0 mul, 2 add
      sink(E4Sol{t.a, t.b, t.c, t.d});
  }
}
```

The norm id is not a parameter: on the host the meet runs inside the DFS, which
already knows which norm it is on, so the sink closes over it. The device sink
will do the same, constructing its `MeetHit` from a `norm_id` it was built with.

Two things this buys beyond tidiness:

- **A8 falls out for free.** Precomputing `sq[]` once per norm turns every key
  from two multiplies + an add into a single add. At $D = 4$ that is 12
  multiplies → 4, in the loop that is 75% of all work.
- `sq[D]` is `D × 8` B **per thread** — 32 B at $D=4$, pure registers. At
  $D = 192$ it is 1.5 KB and would spill, which is precisely why the dispatch
  below sends large $D$ elsewhere.

### 7.1 Bucket dispatch

$D(cq) = 2D(c)$, so the bucket is known *before* generation — sort the cofactor
table by $D$ once and batches are homogeneous with no stream compaction.

```cpp
template <MeetKey Key, HitSink S>
void meet(int D, const norm_t* nu, MeetScratch<Key>& scratch, S&& sink) {
  switch (D) {                       // 97.0% of norms land in the tables
    case 4: return meet_direct<4, Key>(nu, sink);
    case 5: return meet_direct<5, Key>(nu, sink);   // odd D is real: 4225
    case 6: return meet_direct<6, Key>(nu, sink);
    case 7: return meet_direct<7, Key>(nu, sink);
    case 8: return meet_direct<8, Key>(nu, sink);
    default: return meet_hashed<Key>(D, nu, scratch, sink);  // 3%, 74% of work
  }
}
```

`meet_hashed` computes $D(D{-}1)/2$ keys, sorts them, and reads off equal-key
groups — $O(D^2)$ rather than $O(D^4)$, and **no disjointness check on the
candidate path**, because $\nu$ is injective on a canonized point set. (There is
one on the *emission* path, which runs only on an actual key collision and so is
free; it exists so that a caller passing a non-canonized set gets nothing rather
than garbage.)

**Where the $D\le8$ cut comes from.** Swept on the real $B=10^7$ corpus, single
thread: all-hashed 51.4, direct to $D\le4$ 36.6, $\le6$ 34.6, $\le7$ 35.1,
**$\le8$ 33.7**, $\le10$ 34.4 ns/norm. Worth recording that a microbenchmark
repeating one $D$ says the opposite — it has hashed beating `meet_direct<8>` by
27% — because a fixed $D$ makes the sort's branches perfectly predictable and
never resizes the scratch. The mixed corpus is the workload.

**What each idea bought,** measured over every admissible norm with $\nu$ values
precomputed, against two baselines:

| $B$ | oracle: exact, no hoist | A8: hoisted `sq[]`, trunc key | A8+C9: dispatched | hoist | bucket |
| --- | --- | --- | --- | --- | --- |
| $10^6$ | 43.9 ns/norm | 32.0 | 20.4 | 1.37× | 1.57× |
| $10^7$ | 75.0 | 55.3 | 34.2 | 1.36× | 1.62× |
| $10^8$ | 123.9 | 100.3 | 48.2 | 1.24× | 2.08× |

Hoisting is a flat ~1.3×. Bucketing is the one that scales, tracking the $D>8$
work share and sitting near its Amdahl ceiling at each bound. The ~20× projected
at $B = 10^{12}$ is consistent with this trend but is *not* measured.

## 8. The DFS spine

The generator is a state machine over the cofactor tree, carrying point sets
incrementally. Written once; the host walks it directly, the device gives one
block per subtree with the ancestor path in shared memory.

```cpp
class Spine {
public:
  explicit Spine(const FactorBase& base);
  void push(factor_index_t i);      // extend by one prime copy; O(D) multiplies
  void pop();                       // a decrement
  std::span<const Gauss> points() const;
  std::span<const factor_index_t> factors() const;
  norm_t value() const;
private:
  std::array<factor_index_t, kMaxSpineDepth> path_;
  std::array<int, kMaxSpineDepth + 1> size_, offset_;
  std::vector<Gauss> buffer_;       // levels packed end to end, one allocation
};
```

**The relation, stated carefully, because a plausible version of it is wrong.**
$U(cq) = \{\mathrm{canon}(zq), \mathrm{canon}(z\bar q) : z \in U(c)\}$,
*deduplicated*. $U(c)$ is not the set of Gaussian integers of norm $c$ — it is
that set modulo units **and conjugation**. The identity holds because the full
set $W(c)$ is closed under conjugation, so $\overline{W(c)q} = W(c)\bar q$.

The deduplication is not a single collision. At $c = 25$, $q = 2+i$: the parent
$U(25) = \{(5,0),(4,3)\}$ produces $(10,5), (10,5), (10,5), (11,2)$ — a
**three-way** collapse — giving $U(125) = \{(10,5),(11,2)\}$. Anything that
assumes "at most one duplicate, from the self-conjugate point" returns a point
set that is too large, and a point set that is too large is one on which $\nu$
is not injective, which is the precondition the whole meet rests on. So:
canonize, sort, unique, at every level.

**Budgets (corrected).** Exact, from a search over products of split primes
$\le B$, and agreeing with the exhaustive enumeration wherever one exists
($D_{\max} = 16$ at $10^6$, $36$ at $10^8$, matching §2):

| $B$ | max prime copies (depth) | max divisors | max $D$ |
| --- | --- | --- | --- |
| $10^6$ | 8 | 32 | 16 |
| $10^8$ | 11 | 72 | 36 |
| $10^{10}$ | 14 | 160 | 80 |
| $10^{12}$ | **17** | 384 | 192 |

An earlier version of this section budgeted `kMaxDepth = 13` for $10^{12}$. That
is wrong, and the error is instructive: depth is driven by the *smallest* split
prime, not by how many distinct primes fit. $5^{17} = 7.6\times10^{11}$ is an
admissible norm below the design bound with seventeen prime copies — and only
nine points. $D$ is driven by the opposite extreme. Shipped constants are
`kMaxSpineDepth = 20`, `kMaxPointCount = 192`; both are enforced with a throw,
never a truncation.

**Measured** (single thread, `-O2 -march=native`, every admissible norm):

| $B$ | DFS walk only | rebuild per norm | spine | point sets alone | speedup |
| --- | --- | --- | --- | --- | --- |
| $10^6$ | 12.9 ns/norm | 139.7 | 63.0 | 126.7 → 50.1 | 2.22× |
| $10^7$ | 11.5 | 135.4 | 50.8 | 123.9 → 39.3 | 2.67× |
| $10^8$ | 10.8 | 137.0 | 53.7 | 126.2 → 42.9 | 2.55× |

Generation plus meet — the campaign's actual inner loop — is **115 ns/norm at
$B = 10^8$**, against 137 for the reference's point-set rebuild alone.

`NormView` also carries $\nu$, computed at emission rather than in `push` (which
runs on every node of the cofactor tree, most of which are not admissible
norms). Measured honestly, this is a wash for speed: it costs ~5 ns/norm for a
consumer that does not want $\nu$ and saves about the same for one that does.
It is there so the meet reads the walk's own buffer, rather than every caller
re-deriving the point-order convention its injectivity precondition depends on.

## 9. Backends

As shipped, the campaign is not templated on a backend — there is one, the host,
and `std::thread` over an atomic segment counter rather than OpenMP. A device
backend would slot in at the same seam (`Walker` + `meet<Key>` + `certify`), but
inventing the abstraction before the second implementation exists would be
guessing at its shape.

```cpp
class Campaign {
public:
  explicit Campaign(CampaignOptions opts);
  CampaignStats run();              // segmented, checkpointed, resumable
};
```

**The exception, deliberately preserved:** `meet_oracle` and `Enumerator` in
`cpu/enumerate/` stay separate, naive, exact-arithmetic implementations. They
are not the production path with different parameters — they are independently
written oracles, and the cross-check is worthless if the two share code.

### 9.1 Segments, and why the batch size cannot be a constant

The unit of work, scheduling and checkpointing is a **segment**: all norms whose
largest prime factor is one particular prime. Segments partition the space, so
progress is fully described by "segments $0..k$ are done".

Two different balance questions, with opposite answers, and getting them
confused costs either the parallelism or the checkpoint:

- **Across threads, segments are well balanced.** At $B = 10^8$ the largest
  single segment is 2.2% of total cost, against 89.4% if you segment by the
  *smallest* prime instead (the natural ascending DFS order — only 42 non-empty
  segments, and a speedup ceiling of 1.1× on any machine). Dynamic pull
  scheduling over 4.3M segments needs nothing cleverer.
- **Across consecutive indices, segments are wildly unbalanced.** Segment cost
  falls steeply with index, because small primes admit far more smooth
  cofactors. Measured over $B = 10^{10}$ with a fixed batch of 4096:

  | batch | segments | norms | share |
  | --- | --- | --- | --- |
  | 0 | 0–4,096 | 65,199,984 | **39.4%** |
  | 1 | 4,096–8,192 | 15,572,736 | 9.4% |
  | 2 | 8,192–12,288 | 8,562,385 | 5.2% |
  | … | 1,056 batches | … | … |

  So the first checkpoint lands 39% of the way in. Geometric batch growth does
  not fix it — the concentration is *inside* the first block, not spread over
  the first few. The shipped answer is to **adapt the batch size to measured
  wall time**, rescaling after each batch towards a target interval (default
  5 s, clamped to a 4× swing per step). No model of the cost curve, works at any
  bound. Observed ramp at $B = 10^{10}$ with a 2 s target: 64 → 227 → 526 →
  1,482 → 5,219 → 21,690 → 264,744 segments.

### 9.2 The commit order

```text
1. write the batch's rows to the catalog
2. fflush + fsync                          <- rows are durable
3. write checkpoint.tmp, fsync it
4. rename checkpoint.tmp -> checkpoint     <- atomic on POSIX
```

A crash between 2 and 4 leaves rows the checkpoint does not know about, so
resume **truncates** the catalog to the checkpointed byte offset before
appending. Truncating is what makes replay idempotent; appending blindly would
duplicate a batch. Rows are written back in segment order, so the catalog is
byte-identical regardless of thread count or batch boundaries — which is what
lets "resumed" and "uninterrupted" be compared at all.

### 9.3 The B = 10^10 campaign

Measured end to end, 16 threads, factor base loaded from
`data/factorbase_1e12.bin`:

| | |
| --- | --- |
| segments | 4,323,837 |
| admissible norms | 165,632,402 |
| factor base (load, not sieve) | 0.05 s |
| enumeration | **17.0 s** (9.7 M norms/s; 11.0× on 16 threads) |
| meet candidates | 88 |
| rejected, key collision mod $2^{64}$ | 0 |
| rejected, imprimitive | 80 |
| catalogued, all certified | **8** |

Resumability was tested by `SIGKILL` — not a clean stop — eight times in
succession during a live run; the final catalog and checkpoint are byte-identical
to the uninterrupted run. All eight catalog rows were then independently
reproduced in Python, including the two non-squarefree norms
($37^2\cdot101\cdot401$ and $53\cdot73^2\cdot857$), and an independent Python
exhaustive search to $10^6$ agrees on both the norm count (11,083) and the
solution set.

---

# Part IV — Nuances and pitfalls

Freestyle, roughly in order of how much damage each does.

**Silent-miss failures are the only ones that matter.** An exclusion claim that
is wrong produces *nothing*, forever, and looks exactly like a correct one. Every
design choice below that trades safety for speed must be re-examined against
this. No function in the enumeration path may return "empty" to mean "I failed".

**A composite can have a perfectly good decomposition.** $85 = 5\cdot17$,
$125 = 5^3$, $221 = 13\cdot17$ are sums of two squares, so Cornacchia closes and
returns a valid decomposition of a non-prime. Inferring primality from "the
descent worked" is wrong, and these are exactly the composites a caller passes
by accident. Deterministic Miller–Rabin, not cleverness.

**The $B/65$ cutoff is exactly tight.** At $B = 10^6$ the largest prime factor
occurring is 15,377 against a cutoff of 15,384; the extremal norm is
$5\cdot13\cdot15377 = 999{,}505$. Seven units of slack. Any "simplification" of
that bound silently deletes solutions.

**Primitivity: reject only when all four points are imprimitive.** A quadruple
with *some* points divisible by $p$ is genuinely primitive — that is exactly the
case where one pair uses $p = \mathfrak{q}\bar{\mathfrak{q}}$ while another uses
$\mathfrak{q}^2$. The obvious rule ("any imprimitive point ⇒ duplicate") is
wrong and loses real solutions. Imprimitive share is $(b{-}1)/(b{+}1)$.

**Non-squarefree norms need branching over each *copy* of a repeated prime.**
Iterating over distinct primes realizes only exponent 0 or $b$; the full range
$0 \le c \le b$ is what Prop. 5.1 requires. ~20% of admissible norms are
non-squarefree and they carry primitive solutions.

**Canonization is not deduplication — skip it and $\nu$ stops being injective.**
$\nu(\bar z) = \nu(z)$, so on the raw products of point-set construction $z$ and
$\bar z$ carry the same value. The disjointness test is omitted from every meet
*because* $\nu$ is injective, and it is canonization that makes it so; without
it, $\{z_1,z_2\}$ and $\{z_1,\bar z_2\}$ have equal keys and share a point, and
the kernel reports a non-solution. The tempting counter-argument — "the key uses
$\nu^2$, so the sign flip $\nu(iz) = -\nu(z)$ does not matter" — is true and
irrelevant: the danger is conjugation, not the units. See §2.1.

**$\nu$'s injectivity is per-norm.** It licenses dropping the disjointness test
*within* a point set, because $m^2+n^2 = L_1$ and $m^2-n^2 = \nu$ determine
$(m,n)$. It says nothing across norms. Any future kernel that batches points
from different $L_1$ into one hash table reintroduces the need for the check.

**Narrowing `GaussT<int32_t>` does not narrow `norm()` or `nu()`.** They reach
$10^{12}$ while a coordinate reaches $10^6$. The `wide` result type exists for
this; removing it "because coordinates are 32-bit" is a silent overflow in the
one place that is hardest to notice.

**`__int128` has no hardware divide.** Fine for the confirm path (rare), fatal
if it creeps into a loop. Keep it behind `ExactKey` and out of the meet.

**`__launch_bounds__(N, 1)` is not "unconstrained".** It promises one resident
block, so nvcc spends more registers and occupancy *drops* — measured 66 regs /
58% versus 56 regs / 75% with no second argument at all. Always confirm with
`-Xptxas -v` rather than reasoning about it.

**Pinned allocation must be pooled.** `cudaHostAlloc` per call took the demo's
GPU stage from 0.92 ms to ~4 ms. Invisible at demo size, ruinous at campaign
scale. This is currently *unfixed* and needs the session object.

**`atomicAdd` requires exactly `unsigned long long`.** On LP64 `std::uint64_t`
is `unsigned long` — a distinct type. Getting this "right" with a cast
reintroduces the pointer-punning the named types removed.

**`#ifdef __CUDACC__` members are an ODR violation.** Benign here — identical
data members, so layout and by-value kernel arguments are unaffected, and it is
what Thrust and CUB do. The invariant that keeps it benign: **nothing inside the
guard may touch layout.** A static member or a virtual would break it silently.

**`double` runs at 1/64 rate on sm_89.** Every kernel is integer-only. A stray
floating-point `sqrt` in a device path is a 64× cliff.

**Atomic append makes hit order nondeterministic.** Any CPU/GPU comparison must
sort first. A test that compares raw output will flake.

**The exemplar sits above `kExactKeyBound`.** $L_1 = 33{,}065{,}996{,}717 >
3.04\times10^9$, so it exercises the *filter* path, not the exact one. That
makes it the right regression fixture — and means a bug in confirmation would
show up there and nowhere in the small-$B$ tests.

**Modular filtering loses 6 bits, not 0.** Key differences are always divisible
by $2^6\cdot3^2$, so `uint64` gives 58 effective bits, not 64. Ample
($10^{-6}$ false positives over $4\times10^{11}$ comparisons) but the arithmetic
should be done, not assumed, if the modulus ever changes.

**The least quadratic non-residue is not a random variable.** Searching
$n = 2, 3, 4, \dots$ for a non-residue mod $p$ and capping at $k$ trials "because
half of all residues are non-residues, so it fails with probability $2^{-k}$" is
wrong — those are the *smallest* integers, not random samples. Measured over
every split prime in each factor base, the largest least-non-residue is 41 below
$1.54\times10^6$, 53 below $1.54\times10^7$, and **67 at $p = 48{,}473{,}881$**
below $1.54\times10^8$ — growing ~+13/decade, like $\log^2 p$. A 64-trial cap
justified by the $2^{-64}$ argument rejects that prime as composite. Since
primality is established first, the loop is *guaranteed* to terminate; cap it
only as a runaway guard, generously, and make exhaustion loud.

**`std::vector<bool>` of length $B/65$ is 1.9 GB at $B = 10^{12}$** — but that is
a defect in the monolithic implementation, not a property of the problem. A
segmented sieve holds 0.1–0.2 MB resident at *any* bound (measured). Do not
conclude from the 1.9 GB that the sieve is a barrier; see §3.1.

**Store the factor base as $(a,b)$, never $(p,a,b)$.** $p = a^2+b^2$ is
recoverable, so the prime is pure redundancy: 8 B/prime instead of 16, and the
whole $B = 10^{12}$ base is 2.7 GB rather than 5.2. It is an `mmap`'d artifact
streamed sequentially — it never needs to be resident on either host or device.

**`Norm::primes` aliases live DFS state.** Deliberate — a value type is a heap
allocation per norm in a $2\times10^{10}$-iteration loop — but any callback that
stores the reference gets garbage. Copy in the callback.

**Point-set construction is two allocations per norm.** At 154 ns/norm this is
already the dominant CPU cost, and the allocations are most of it. The spine
(§8) removes them by construction; until then, do not micro-optimize the meet.

**Checkpoint at segment boundaries, and pre-register $B$.** A campaign that
cannot resume will not finish, and a bound chosen after seeing partial results
is not an exclusion claim.

**Two independent implementations, and keep them independent.** The value of the
CPU/GPU cross-check is that neither can hide the other's bug. Sharing `meet<P>`
between backends is right; sharing it with `meet_reference` would quietly turn
the fixture into a tautology.

**The family self-test is the one that finds coverage holes.** Re-finding
Bremner's two families *including the non-squarefree $L_1(k)$* is what exposes an
A3/A3.1 gap. Norm counts and planted collisions will not.
