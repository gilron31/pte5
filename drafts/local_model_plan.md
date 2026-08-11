# G1 — The Local Model for $\sigma_p$: research plan and implementation plan

**Status of this document.** A plan, not a result. Every number below marked
**[measured]** comes from a run recorded here; everything else is a proposal or
a target. The one thing it does claim is a *diagnosis*: the two models
`gpu_implementation_plan.md` §2.2 rejects are not two competing answers, they
are the same wrong question asked twice, and naming the right question makes all
three of §6.6's deliverables fall out of one computation.

---

## 0. The question, and why the two obvious answers both fail

`paper.md` §6.6 wants, for each prime $p$, a local density $\sigma_p$ feeding the
singular series of §8, plus sieve wheels, plus a possible nonexistence
certificate. The obvious computation — count $(L_1,\dots,L_d) \in \mathbb{F}_p^d$
whose tower splits — has two variants and neither is $\sigma_p$:

| model | what it counts | why it is wrong |
| --- | --- | --- |
| **distinct roots** | towers whose $2^d$ roots are distinct mod $p$ | **Proves a false theorem.** Zero at 25 of the 29 primes $\le 113$, yet a *known* $E_4$ exists. Its 16 roots are distinct over $\mathbb{Z}$ and collide mod $p$ at 16 of those primes |
| **with multiplicity** | towers splitting mod $p$, repeats allowed | **Counts only degeneracies at small $p$.** Once $A$ collapses to $\{0\}$ the sole completion is $L_d = 0$, i.e. $E_d = E_{d-1}^2$ — a tower that is not a solution |

> **The diagnosis.** Both are the $k = 1$ truncation of a $p$-adic quantity, and
> the truncation destroys exactly the information that distinguishes them. A
> genuine solution is a $\mathbb{Z}$-point, hence a $\mathbb{Z}_p$-point, of the
> **nondegenerate** locus. Nondegeneracy is not a mod-$p$ condition — it is the
> statement that a node value is nonzero *in $\mathbb{Z}_p$*, i.e. that its
> valuation is **finite**. Mod $p$ you can see "zero" and "nonzero" but not
> "finite"; that distinction first becomes visible mod $p^k$, and only for $k$
> larger than the valuation actually occurring.

Everything below follows from taking that seriously.

---

## 1. The proposed object

**Towers over a ring.** For a commutative ring $R$, a *depth-$d$ tower* over $R$
is $(L_1,\dots,L_d) \in R^d$ together with a value at every node of the binary
tree of nested radicals: $u^2 = L_d$ at the top, and $w_\pm^2 = L_j \pm v$ for
each node value $v$ at level $j{+}1$, down to the roots $r^2 = L_1 \pm(\cdot)$.
This is `paper.md` §2's condition verbatim, and it is well defined over
$\mathbb{Z}/p^k$ without needing "splits into linear factors" to mean anything
over a non-domain. It is also exactly what the §2.1 recursion computes:

$$T_d = \{\,t^2 - L_d : t \in A_{d-1}\,\},\qquad A_j = \{t \in T_j : -t \in T_j\},\qquad T_0 = R,$$

with the tower splitting iff $0 \in T_d$. **The recursion is unchanged mod
$p^k$** — only the ring changes.

**Degeneracy.** A tower is degenerate iff some node value is $0$: a zero node
forces its two children to coincide, collapsing the subtree and producing
repeated roots. Over $\mathbb{Z}$, "nondegenerate" is `paper.md`'s nontriviality.
Over $\mathbb{Z}_p$ it is "every node value has finite valuation".

**The valuation profile.** For a tower $\tau$ write
$$\mathrm{ht}_p(\tau) = \max_{\text{nodes}} v_p(\text{node value}) \in \{0,1,\dots,\infty\},$$
and for each precision $k$ define the **valuation floor**
$$V_p(k,d) = \min\{\,\mathrm{ht}_p(\tau) : \tau \text{ a depth-}d \text{ tower over } \mathbb{Z}/p^k\,\},$$
valuations capped at $k$ (indistinguishable from $\infty$ at that precision).

This single quantity carries all three deliverables:

| §6.6 wants | read off as |
| --- | --- |
| **nonexistence certificate** | $V_p(k,d) \to \infty$ in $k$ ⟹ every $\mathbb{Z}_p$-tower is degenerate ⟹ **no $E_d$ over $\mathbb{Z}$** |
| **forced divisors** | $V_p(k,d) \to v_0 < \infty$ ⟹ every $E_d$ has a node value divisible by $p^{v_0}$ — the §4 divisibility ladder, *computed* rather than proved |
| **local density $\sigma_p$** | $\lim_k p^{-k\dim} N_k$ restricted to $\mathrm{ht}_p \le v$, summed over $v$ |

The compactness argument for row 1 is the only rigorous inference here and it is
short: a nondegenerate $\mathbb{Z}_p$-tower has some finite $\mathrm{ht}_p = v_0$,
and its reduction mod $p^k$ is a tower with $\mathrm{ht}_p \le v_0$ for every
$k > v_0$; so if for every $v$ there is a $k$ with $V_p(k,d) > v$, no such tower
exists. **This is the branch that can produce a theorem.**

**Density, and an honest limit.** $\sigma_p$ is the $p$-adic density
$\lim_k p^{-k(N-c)} N_k$ in the variables of `paper.md` §8.1 — which are the
*root* coordinates, not the $L$'s. §6.6 currently counts $L$-tuples; the two
differ by a fibre count that is not constant across the degenerate locus, so
**fixing the coordinates is a task, not a detail** (§4, task T3). By Igusa's
theorem the Poincaré series $\sum_k N_k p^{-Nk} t^k$ is *rational*, so the
sequence satisfies a linear recurrence and the limit can be extracted exactly
from finitely many terms rather than guessed from a plateau.

> **What the singular series cannot do here.** $V_5 \subset \mathbb{P}^{15}$ has
> 16 variables and total degree 34; the circle method needs the opposite
> inequality by a wide margin, and Birch's condition on the singular locus fails
> too. $\prod_p \sigma_p$ is therefore a *heuristic weight* for §8 and can never
> be upgraded to a theorem by this route. **Only the obstruction branch is
> capable of proving something**, which is an argument for computing it first
> and separately.

---

## 2. First measurements

A prototype of the mod-$p^k$ recursion with a cap on $\mathrm{ht}_p$
(`drafts/` scratch, to be productionised as T1). Entries are $V_p(k,d)$; a value
equal to $k$ means "nothing nondegenerate is visible at this precision".

**[measured]**

| $p$ | $k$ | $d{=}3$ | $d{=}4$ | $d{=}5$ | towers mod $p^k$ (uncapped) |
| --- | --- | --- | --- | --- | --- |
| 3 | 1 | 1 | 1 | 1 | 3 |
| 3 | 2 | 2 | 2 | 2 | 7 |
| 3 | 3 | **2** | 3 | 3 | 57 |
| 3 | 4 | **2** | 4 | 4 | 385 |
| 3 | 5 | **2** | 5 | 5 | 3,423 |
| 5 | 2 | **1** | 2 | 2 | 51 |
| 5 | 3 | **1** | **2** | 3 | 3,755 |
| 7 | 2 | **1** | 2 | 2 | 211 |

Depth 3 stabilises immediately; depth 4 stabilises at $p = 5$; depth 5 tracks
$k$ everywhere so far. **That shape is what a nonexistence proof would look
like — and it is also what "not deep enough yet" looks like**, which is why the
next measurement is the control, not the extension.

**The control, and it validates the prototype.** The §7 exemplar $E_4$ is a real
depth-4 tower, so its reduction must appear in the table once $k$ is large
enough. Its node valuations, computed directly: **[measured]**

| level | $v_2$ | $v_3$ | $v_5$ | $v_7$ |
| --- | --- | --- | --- | --- |
| 1 (roots) | 0 | 0, 1 | 0, 1 | 0 |
| 2 | 2 | 0 | 0 | 0 |
| 3 | 6 | 3, 4 | 3, 4 | 1 |
| 4 ($\sqrt{L_4}$) | **14** | **6** | **6** | **2** |

So $\mathrm{ht}_3 = 6$ for a *known* solution — and the prototype's failure to
find a nondegenerate depth-4 tower mod $3^5$ is therefore **correct, not a bug**:
$5 < 6$. It also reproduces §4 independently ($v_3(C'_8/2) = 3 \le 6$,
$v_2(C'_8/2) = 5 \le 14$, $v_3(C'_4/2) = 2 \le 3$), and the 2-adic column shows
the ladder as a clean progression $0, 2, 6, 14 = 2^{\ell}-2$ that `paper.md` §2's
normalizing constants $2^{2^k-1}$ predict.

> **This fixes the required precision, which was the missing number.** You cannot
> conclude anything about depth $d$ at prime $p$ until $k$ exceeds the valuation
> a real solution needs. §4's $C'_{16} = 2^{11}3^6 5^4 7^3 11^2 13^2 \cdot 17
> \cdots 53$ gives the floor for depth 5 directly: $k \ge 12$ at $p=2$, $k \ge 7$
> at $p=3$, $k \ge 5$ at $p=5$, $k \ge 4$ at $p=7$, $k \ge 3$ at $p = 11, 13$,
> $k \ge 2$ at $p \le 53$. **Every published depth-5 local claim computed at
> $k=1$ — including the plan's own table — is below this floor.**

---

## 3. Calibration ladder

Three independent known answers, in increasing order of strength. Nothing about
depth 5 is trustworthy until all three pass.

1. **Depth 3 against [TLT].** `paper.md` §8.2 quotes the *rigorous* count
   $N_{\mathrm{sym}}(H) = \frac{4\log 2}{3\pi^2}H^3\log H + O(H^3)$. Assembling
   $\prod_p \sigma_p$ at $d=3$ with the archimedean factor must reproduce that
   constant. This is the only place in the project where a density heuristic can
   be checked against a theorem, and it tests the *whole pipeline* — coordinates,
   normalization, convergence — not just the local factors.
2. **Depth 4 against the catalog.** 8 certified solutions at $B = 10^{10}$, exact
   admissible-norm counts at $10^6/10^7/10^8$, and Bremner's two infinite
   families. The families live on a special subvariety and are *expected* to
   break a naive density prediction (§8.1: the generic count is $H^{-2}$ yet
   families exist) — so the test is that the model's error is the *family
   contribution*, identified as such, not an unexplained factor.
3. **Forced divisors against [CMSV].** $V_p(k,4)$ must stabilise at a value
   consistent with $C'_8 = 2^6 3^3 5^2 7^2 \cdot 11 \cdot 13$, and $V_p(k,5)$
   with $C'_{16}$. This is the sharpest test because it is *exact arithmetic
   against exact arithmetic*, with no heuristic in between — and if the model
   produces a divisor **larger** than CMSV's, that is a new theorem about PTE
   constants, obtained by machine.

---

## 4. Work items

| | task | output | depends |
| --- | --- | --- | --- |
| **T1** | Mod-$p^k$ tower DP with valuation cap, in C++, on $p^k$-bit masks | $V_p(k,d)$, $N_k$, state counts | — |
| **T2** | Brute-force fixture: enumerate all $(L,\text{tree})$ mod $p^k$ for tiny $(p,k,d)$ and compare | the oracle §6 of the GPU plan demands | T1 |
| **T3** | Fix the coordinates: define $\sigma_p$ in root coordinates, derive the $L$-count as a consequence, state the fibre correction | an unambiguous $\sigma_p$ | — |
| **T4** | Exemplar and family regression: every known $E_4$ must appear as a tower mod $p^k$ with its measured valuation profile | the control that makes a null result meaningful | T1 |
| **T5** | Poincaré series: fit a linear recurrence to $N_k$, verify, extract $\lim_k$ exactly | $\sigma_p$ as a rational number | T1, T3 |
| **T6** | Push depth 5 to the §4 floor: $k \ge 7$ at $p=3$, $\ge 5$ at $p=5$, $\ge 12$ at $p=2$ | the obstruction test, honestly conducted | T1, T4 |
| **T7** | Assemble $\prod_p \sigma_p$; run calibration 1 and 2 | the §8 constant, or a refutation | T3, T5 |
| **T8** | Sieve wheels: residues of $L_1$ mod $p^k$ carrying no nondegenerate tower | a pre-sieve for S-A / S-C | T1 |

**Order.** T1 → T2 → T4 → T6 is the obstruction branch and is the only one that
can end in a theorem; it should be finished before T3/T5/T7, which are the
density branch and can only ever end in a heuristic.

---

## 5. Implementation

**Where.** `code/native/src/cpu/tower/` — a new CPU app directory, mirroring
`champion/`: `tower.h` (the DP over a modulus), `tower.cc`, `tower_test.cc`.
Pure CPU, no GMP needed (everything is mod $p^k \le 2^{15}$), so it links
`pte_core` alone.

**The state, and why this is tractable.** The DP state is the set
$T \subseteq \mathbb{Z}/p^k$ as a $p^k$-bit mask. The transition
$T' = \{s - L : s \in S\}$, $S = \{t^2 : t \in A\}$, is a **cyclic shift of the
mask by $L$** — the same structure as `gpu_implementation_plan.md` §2.3
describes mod $p$, and it survives verbatim mod $p^k$. So one level costs
$|\text{states}| \times p^k \times p^k/64$ word-operations.

Sizing from the measured state counts (states per level, $p=3$): $[1,3,3,4,4]$
at $k{=}1$, $[1,27,55,28,28]$ at $k{=}3$, $[1,243,730,244,244]$ at $k{=}5$ —
**states grow like $O(p^k)$, not $O(p^{2k})$**, which is the fact that makes the
whole thing feasible. At the §4 floor:

| $p^k$ | mask | states (extrapolated) | state table | word-ops/level |
| --- | --- | --- | --- | --- |
| $3^7 = 2{,}187$ | 274 B | ~6,600 | 1.8 MB | $2.4\times10^{8}$ |
| $5^5 = 3{,}125$ | 391 B | ~9,400 | 3.7 MB | $4.6\times10^{8}$ |
| $2^{12} = 4{,}096$ | 512 B | ~12,000 | 6.1 MB | $7.9\times10^{8}$ |
| $5^7 = 78{,}125$ | 9.8 KB | ~235,000 | 2.3 GB | $2.3\times10^{13}$ |

The first three rows are **seconds** on one core. The last is the honest wall:
matching the exemplar's $v_5 = 6$ needs $5^7$, and that is a GPU job or a
better algorithm. Measure the state growth at each step rather than trusting the
extrapolation — it is the one quantity that could kill this.

**Dedup.** Distinct states must be recognised; hash the mask (64-bit
fingerprint, verify on collision) into an open-addressing table — the same
structure `gpu/champion` already has, and the same one §2.3 of the plan
specifies for the device version.

**Tests.** (a) T2's brute force for $p^k \le 5^2$, $d \le 3$, every prefix, exact
set equality — this is the fixture the GPU plan's §6 says must ship. (b) The
$k=1$ results must reproduce `drafts/verify_gpu_plan.py`'s table exactly
(13/26/14/1 states at $p=13$, the depth-5 counts 0/78/11/207/371 at
79/89/97/113). (c) T4's exemplar profile, as a table of expected valuations.
(d) $V_p(k,d)$ monotone non-decreasing in $k$ and in $d$ — a cheap invariant
that catches most cap-handling bugs.

**The device question, which the correction changes.** `gpu_implementation_plan.md`
§1 ranks S-D behind S-E partly because P5 measured the bitset rotate at
1.5e10 word-op/s at $p = 1024$ against 1.1e11 at $p = 65536$ — "the device is at
its worst exactly where S-D lives". **Under the corrected model S-D no longer
lives there**: the modulus is $p^k \approx 2\times10^3$ to $8\times10^4$, which
is the regime where P5 measured the device at its *best*. G5's conditional
should be re-read after T1.

---

## 6. Prior art

The problem sits at the intersection of three literatures, and the useful thing
is that each supplies a different half.

- **$p$-adic densities and point counts mod $p^k$.** The generating function
  $P(t) = \sum_k N_k p^{-Nk}t^k$ is Igusa's Poincaré series, and its
  **rationality** (conjectured by Borevich–Shafarevich, proved by Igusa via
  resolution of singularities) is what turns "does the density stabilise?" into
  "fit a linear recurrence and verify". This is the tool for T5, and it also
  says the $V_p(k,d)$ growth in §2 should be *eventually linear or eventually
  constant* — a dichotomy that is precisely the theorem/forced-divisor split.
  ([survey](https://eudml.org/doc/110143),
  [computational side](https://msp.org/obs/2020/4-1/obs-v4-n1-p13-s.pdf))
- **Arboreal Galois representations.** An $E_d$ tower is a totally-split
  preimage tree for the non-autonomous quadratic composition $y^2 - L_i$ —
  `paper.md` §9 says so explicitly. The arboreal literature (Jones's survey;
  Boston–Jones on settled elements) studies exactly the density of primes for
  which such trees split, via Chebotarev on the tree's splitting field, and
  supplies the *other* cross-check: for a **fixed** integer tower the density of
  good primes is $1/|G|$ for an explicitly computable $G$. Our exemplar's
  16-roots-collapse-at-16-of-29-primes observation is a coarse instance; making
  it exact is an independent test of T4.
  ([survey](https://www.researchgate.net/publication/260366988_Galois_representations_from_pre-image_trees_an_arboreal_survey),
  [uniformization](https://arxiv.org/pdf/1111.3607))
- **PTE required divisors.** Kleiman, Rees–Smyth, Caley, Filaseta–Markovich and
  [CMSV] compute $C_n$, $C'_n$ by hand, "by local and $2$-adic arguments"
  (`paper.md` §4). **That is this computation, done by hand, for the top node
  only.** T6 is its mechanisation and extension to every node of the tree, and
  calibration 3 is the check that the machine agrees with the humans.
  ([CMSV](https://arxiv.org/abs/2304.11254),
  [Filaseta–Markovich via the review](https://www.mdpi.com/2227-7390/7/3/227),
  [TLT](https://arxiv.org/pdf/2606.07735))

What none of them supplies, as far as these searches went, is the specific
object: the $p$-adic density of *non-autonomous* totally-split quadratic towers
with the degenerate locus excluded. If T6 produces a divergent $V_p(k,5)$, that
is a new result and the write-up should say what it would take to make it a
proof — namely the eventual-periodicity of the DP, which is where Igusa
rationality re-enters.

---

## 7. Stopping rules

Written in advance, because this is the branch where wishful reading is
cheapest:

- **A depth-5 null at $k$ below the §4 floor means nothing** and must not be
  reported. The floor is $k \ge 12, 7, 5, 4, 3, 3, 2$ at $p = 2, 3, 5, 7, 11,
  13, \le 53$.
- **A depth-4 null at any $k$ is a bug**, since solutions exist — T4 is the
  tripwire.
- $V_p(k,5)$ growing 1:1 with $k$ **past** the floor, at a prime where
  $V_p(k,4)$ has stabilised, is the signal. Below that it is indistinguishable
  from insufficient precision.
- The density branch is a heuristic and stays labelled one, whatever it produces
  (§1's note on the circle method).
