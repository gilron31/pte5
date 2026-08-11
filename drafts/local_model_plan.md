# G1 — The Local Model for $\sigma_p$: research plan and implementation plan

**Status of this document.** A plan with its first three items built. Every
number marked **[measured]** comes from a run recorded here; everything else is
a proposal or a target. The claim it opens with is a *diagnosis*: the two models
`gpu_implementation_plan.md` §2.2 rejects are not two competing answers, they
are the same wrong question asked twice, and naming the right question makes all
three of §6.6's deliverables fall out of one computation.

**Built (T1, T2, T4).** `code/native/src/cpu/tower/` — the mod-$p^k$ recursion
with the valuation cap (`tower.h`/`tower.cc`), the tree-condition oracle it is
checked against (`tower_oracle.h`/`.cc`), and `tower_test.cc`; the driver is
`apps/tower_demo.cc`. 10/10 tests in 51 ms, including exact agreement with
brute force at every cap for ten $(p,k,d)$ combinations. §2 below is now
measurement rather than proposal.

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

## 2. Measurements

All from `apps/tower_demo`, single core. Entries are $V_p(k,d)$; a starred value
equals $k$ and means "nothing nondegenerate is visible at this precision" — a
non-answer, not a null. An unstarred value is stable and carries the inference
of §2.1.

**[measured]**

| $p^k$ | $d{=}2$ | $d{=}3$ | $d{=}4$ | $d{=}5$ | states/level ($d=5$, uncapped) | s |
| --- | --- | --- | --- | --- | --- | --- |
| $2^{6}$ | 2 | **5** | 5\* | 5\* | 1 64 65 65 65 | 0.00 |
| $2^{10}$ | 2 | **5** | 9\* | 9\* | 1 1024 3713 1025 1025 | 0.49 |
| $2^{12}$ | 2 | **5** | 11\* | 11\* | 1 4096 19073 4097 4097 | 33.7 |
| $3^{5}$ | 1 | **2** | 5\* | 5\* | 1 243 730 244 244 | 0.00 |
| $3^{7}$ | 1 | **2** | **6** | 7\* | 1 2187 8749 4375 2188 | 2.6 |
| $3^{8}$ | 1 | **2** | **6** | 8\* | 1 6561 32806 19684 6562 | 63.2 |
| $5^{4}$ | 0 | **1** | **2** | 4\* | 1 625 1255 1251 626 | 0.06 |
| $5^{5}$ | 0 | **1** | **2** | **4** | 1 3125 9380 9376 6251 | 13.6 |
| $7^{3}$ | 0 | **1** | **2** | 3\* | 1 343 50772 2059 344 | 0.12 |
| $11^{3}$ | 0 | **1** | **2** | 3\* | 1 1331 14642 3994 1332 | 1.5 |
| $13^{3}$ | 0 | 0 | **1** | **2** | 1 2197 8801 4408 4395 | 2.2 |

The exemplar control, now in C++ and reading only the *reduced* tuple — it
recovers $\mathrm{ht}_p$ with no big-integer arithmetic, and agrees with the
direct computation over $\mathbb{Z}$ where the two overlap: **[measured]**

| | $p=2$ | $p=3$ | $p=5$ | $p=7$ | $p=11$ | $p=13$ |
| --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{ht}_p(E_4^{\S 7})$ | 14 | **6** | 6 | **2** | **2** | **1** |
| least $k$ that reveals it | 15 | 7 | 7 | 3 | 3 | 2 |

So the failure to find a nondegenerate depth-4 tower mod $3^5$ is **arithmetic,
not a bug**: $5 < 6$. The 2-adic column is the ladder as a clean $0, 2, 6,
14 = 2^{\ell}-2$, which `paper.md` §2's normalizing constants $2^{2^k-1}$ predict.

### 2.1 What the stable entries prove

If $V_p(k,d) < k$ then **every $E_d$ over $\mathbb{Z}$ has a node value divisible
by $p^{V_p(k,d)}$** — unconditionally, with no hypothesis on the solution. The
argument is three lines: a real tower reduces to a tower mod $p^k$; reduction
sends a node of valuation $v$ to one of capped valuation $\min(v,k)$; so
$\min(\mathrm{ht}_p, k) \ge V_p(k,d)$, and $V < k$ forces $\mathrm{ht}_p \ge V$.

That makes the stable entries directly comparable with `paper.md` Theorem 4.1,
which cascades [CMSV]'s $C'_n$ through the tower and gives
$C'_{2^{d-1}}/2 \mid \sqrt{L_d}$ at the top node — hence
$\mathrm{ht}_p \ge v_p(C'_{2^{d-1}}/2)$. **This is calibration 3, and it is the
sharpest test in the plan: exact arithmetic against exact arithmetic.**
**[measured]**

| | $p=2$ | $p=3$ | $p=5$ | $p=7$ | $p=11$ | $p=13$ |
| --- | --- | --- | --- | --- | --- | --- |
| $d=3$: Thm 4.1 via $C'_4/2 = 18$ | 1 | 2 | 0 | 0 | 0 | 0 |
| $d=3$: computed $V_p(\cdot,3)$ | **5** | 2 | **1** | **1** | **1** | 0 |
| $d=4$: Thm 4.1 via $C'_8/2$ | 5 | 3 | 2 | 2 | 1 | 1 |
| $d=4$: computed $V_p(\cdot,4)$ | ≥11 | **6** | 2 | 2 | **2** | 1 |
| $d=5$: Thm 4.1 via $C'_{16}/2$ | 10 | 6 | 4 | 3 | 2 | 2 |
| $d=5$: computed $V_p(\cdot,5)$ | ≥12 | ≥8 | **4** | ≥3 | ≥3 | **2** |

> **The machine agrees with the theorem everywhere and beats it in four places.**
> Every computed entry is $\ge$ the corresponding required-divisor exponent, which
> is the consistency check passing. Where it is strictly larger — $p=2$ at depth
> 3 (5 vs 1), $p=3$ at depth 4 (**6 vs 3**), $p=5,7,11$ at depth 3, $p=11$ at
> depth 4 — it is a **new forced divisor obtained by computation**, and at $p=3$,
> $d=4$ it is *sharp*: the §7 exemplar attains $\mathrm{ht}_3 = 6$ exactly, so
> $3^6$ cannot be improved.

**Two depth-5 entries are already stable, and both say "no obstruction here".**
$V_5(5,5) = 4$ and $V_{13}(3,5) = 2$ are unstarred, so a depth-5 tower with
finite node valuations exists mod $5^5$ and mod $13^3$ — and both land *exactly*
on $v_p(C'_{16}/2)$. At those two primes the local model neither obstructs an
$E_5$ nor improves on [CMSV]. The obstruction question lives at $p = 2$ and
$p = 3$, which is also where §2.2 says the precision is hardest to reach.

This is `paper.md` §11's "sharper required divisors" item, which asks for a
constant strictly containing $C'_{16}$ and notes that a genuine improvement
"must exploit the *towered* structure specifically". The DP is nothing but the
towered structure. Two caveats keep the claim honest:

- $V_p$ bounds $\max_{\text{nodes}} v_p$, whereas Theorem 4.1 localises to the
  *top* node. The quantities differ, and only the top-node version would tighten
  §4 as stated. Computing it is a small variant — minimise $v_p(t)$ over
  $t \in A_{d-1}$ instead of capping every level — but the inference is **not**
  the same, for the reason in §2.3, and that is what T6b has to solve.
- $d=5$ rows are lower bounds from starred (non-)answers, not results.

### 2.3 Why the *max* is the honest quantity and the *top node* is not, yet

$E_n$ solutions scale: if $(L_1,\dots,L_n)$ is one then so is
$(t^2L_1, t^4L_2, \dots, t^{2^n}L_n)$, because $E_j$ scales by $t^{2^j}$ under
$x \mapsto tx$. Taking $t = p^N$ shows **$\mathrm{ht}_p$ is unbounded above**.

Two consequences, in opposite directions:

- **$V_p$ survives it untouched.** $V_p$ is a *minimum* over towers, and scaling
  only produces towers with larger valuations, so it cannot lower the min. The
  §2.1 inference stands unconditionally. It also explains the sharpness at
  $p=3$: the §7 exemplar is $p$-minimal there, so it sits exactly on the floor.
- **The top-node variant does not.** Its natural statement — "every $E_d$ has
  $v_p(\sqrt{L_d}) \ge F$" — needs the *other* nodes capped, and a real solution
  is only covered once $k$ exceeds its $\mathrm{ht}_p$. With $\mathrm{ht}_p$
  unbounded there is no single $k$ that covers every solution, so what falls out
  is conditional: *every $E_d$ whose node valuations are all $< k$ satisfies
  $v_p(\sqrt{L_d}) \ge F(k)$.*

**The fix is primitivity, and the project already has the notion**
(`src/cpu/enumerate/`). Restricted to primitive towers, $\mathrm{ht}_p$ should be
bounded — and *that* bound, not the DP, is the missing lemma. T6b is therefore
two things: a small code change and a small proof, and the proof is the harder
half. It should not be attempted as an implementation task alone.

### 2.2 The precision floor, revised upward by the data

Theorem 4.1 fixes the *minimum* precision at which a depth-5 claim can mean
anything: $v_p(\sqrt{L_5}) \ge v_p(C'_{16}/2)$, so $k \ge 11, 7, 5, 4, 3, 3$ at
$p = 2, 3, 5, 7, 11, 13$ and $k \ge 2$ at $p \le 53$. **Every depth-5 local claim
computed at $k=1$ — including this plan's own first table — is below that floor.**

But the depth-4 row now says the floor is *optimistic*, and by a measured factor:

| | $p=2$ | $p=3$ | $p=5$ | $p=7$ | $p=11$ | $p=13$ |
| --- | --- | --- | --- | --- | --- | --- |
| $d=4$ floor from $C'_8/2$ | 5 | 3 | 2 | 2 | 1 | 1 |
| $d=4$ truth ($\mathrm{ht}_p$ of the exemplar) | 14 | 6 | 6 | 2 | 2 | 1 |
| ratio | 2.8 | 2.0 | 3.0 | 1.0 | 2.0 | 1.0 |

At the two primes that matter most the required precision is **2–3× the §4
floor**. Extrapolating the same ratio to depth 5 puts the honest target at
$k \approx 15$–$20$ at $p=3$ and $k \approx 25$–$30$ at $p=2$ — that is
$3^{15} = 1.4\times10^7$ and $2^{25} = 3.4\times10^7$, three to four orders of
magnitude past what §5's cost model reaches. **T6 is not a matter of running the
existing code for longer.** §5.1 is how it becomes reachable.

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

| | task | output | depends | state |
| --- | --- | --- | --- | --- |
| **T1** | Mod-$p^k$ tower DP with valuation cap, in C++, on $p^k$-bit masks | $V_p(k,d)$, $N_k$, state counts | — | **done** |
| **T2** | Brute-force fixture: enumerate all $(L,\text{tree})$ mod $p^k$ for tiny $(p,k,d)$ and compare | the oracle §6 of the GPU plan demands | T1 | **done** |
| **T3** | Fix the coordinates: define $\sigma_p$ in root coordinates, derive the $L$-count as a consequence, state the fibre correction | an unambiguous $\sigma_p$ | — | open |
| **T4** | Exemplar and family regression: every known $E_4$ must appear as a tower mod $p^k$ with its measured valuation profile | the control that makes a null result meaningful | T1 | **done** |
| **T5** | Poincaré series: fit a linear recurrence to $N_k$, verify, extract $\lim_k$ exactly | $\sigma_p$ as a rational number | T1, T3 | open |
| **T6** | Push depth 5 past the *revised* floor of §2.2, not the §4 one | the obstruction test, honestly conducted | T1b | blocked |
| **T6b** | Top-node variant: minimise $v_p(\sqrt{L_d})$ instead of capping every level, **plus** the primitivity lemma bounding $\mathrm{ht}_p$ (§2.3) | directly comparable with Thm 4.1, so it can *tighten* $C'_{16}$ | T1 | open — the proof is the hard half |
| **T1b** | Orbit-representative rewrite (§5.1): $O(\text{orbits}\cdot m^2)$ in place of $O(\text{states}\cdot m^2/64)$ | the $m/64$ factor T6 needs | T1 | next |
| **T7** | Assemble $\prod_p \sigma_p$; run calibration 1 and 2 | the §8 constant, or a refutation | T3, T5 | open |
| **T8** | Sieve wheels: residues of $L_1$ mod $p^k$ carrying no nondegenerate tower | a pre-sieve for S-A / S-C | T1 | open |

**Order.** T1 → T2 → T4 → T6 is the obstruction branch and is the only one that
can end in a theorem; it should be finished before T3/T5/T7, which are the
density branch and can only ever end in a heuristic. T1b and T6b are what the
first measurements added: T6b because §2.1 landed next to Theorem 4.1 but not
quite on it, and T1b because §2.2 showed T6 is three orders of magnitude out of
reach without it.

---

## 5. Implementation

**Where.** `code/native/src/cpu/tower/` — a CPU app directory mirroring
`champion/`: `tower.h`/`tower.cc` (the DP over a modulus), `tower_oracle.h`/`.cc`
(the tree condition, per tuple), `tower_test.cc`. It links **nothing**: every
value is a residue mod $p^k \le 2^{20}$ and every state is a bitmask, so no GMP,
no Boost, not even `pte_core`. Worth preserving, because the same DP is the
candidate S-D device kernel and a dependency acquired here would be paid for
again on the device.

**The state.** The DP state is the set $T \subseteq \mathbb{Z}/p^k$ as a
$p^k$-bit mask. The transition $T' = \{s - L : s \in S\}$, $S = \{t^2 : t \in A\}$,
is a **cyclic rotation of the mask by $L$** — the same structure
`gpu_implementation_plan.md` §2.3 describes mod $p$, surviving verbatim mod
$p^k$. Distinct states are interned by 64-bit fingerprint with a full `memcmp`
on a hit, the same shape as `gpu/champion`'s table. One level costs
$|\text{states}| \times p^k \times p^k/64$ word-operations and
$|\text{states}| \times p^k/8$ bytes.

**Measured, and the extrapolation was right.** States per level at $d=5$:
**[measured]**

| $p^k$ | states/level | states / $p^k$ | mask | wall |
| --- | --- | --- | --- | --- |
| $3^{5} = 243$ | 1 243 730 244 244 | 3.0 | 31 B | 0.00 s |
| $3^{7} = 2{,}187$ | 1 2187 8749 4375 2188 | 4.0 | 274 B | 2.6 s |
| $3^{8} = 6{,}561$ | 1 6561 32806 19684 6562 | 5.0 | 821 B | 63 s |
| $2^{12} = 4{,}096$ | 1 4096 19073 4097 4097 | 4.7 | 512 B | 34 s |
| $7^{3} = 343$ | 1 343 50772 2059 344 | 148 | 43 B | 0.12 s |

**States grow like $O(p^k)$, not $O(p^{2k})$** — the fact the whole thing rests
on. But the wall-clock column is the $m^3/64$ law asserting itself: $3^7
\to 3^8$ is $3\times$ the modulus and $24\times$ the time. At that scaling
$3^{15}$, which §2.2 says is the honest depth-5 target at $p=3$, is
$10^{7}\times$ $3^8$'s 63 s. **Two years of core-time.** Hence §5.1.

### 5.1 The orbit collapse, measured

Every state at level $\ge 1$ has the form $S - L$, so each level's state set is
a **union of complete translation orbits**. `tower_demo`'s third table measures
it directly — `misses` counts rotations of a state that are not themselves
states, and it is **0 at every level of every run so far**: **[measured]**

| $p^k$ | states/level | orbits/level |
| --- | --- | --- |
| $3^{4} = 81$ | 1 81 244 82 82 | 1 1 4 2 2 |
| $5^{3} = 125$ | 1 125 255 251 126 | 1 1 3 3 2 |
| $7^{3} = 343$ | 1 343 **50772** 2059 344 | 1 1 **150** 7 2 |
| $2^{8} = 256$ | 1 256 641 257 257 | 1 1 4 2 2 |

The orbit count is usually tiny where the state count is $\Theta(p^k)$ — but
**not universally**, and the exception is worth recording rather than smoothing
over: at $p=7$ the level-2 orbit count grows with $k$ (4, 24, **150** at
$k = 1,2,3$, i.e. $\approx m/2$), while at $p = 3, 5, 11, 13$ it is flat in $k$
(3, 3, 11, 4). $7^4$ is the first row that exceeds the shipped `max_states` and
it does so for this reason.

The speedup does not depend on that, which is the point:

- **now**: for each of $\approx c\,m$ states, materialise all $m$ rotations and
  intern them — $c\,m \cdot m \cdot m/64 = c\,m^3/64$.
- **with orbits**: for each of $\approx c\,m$ states, compute $S$ and canonicalise
  it to (representative, offset) in $O(m)$ by Booth's least-rotation — $c\,m^2$.

The saving is **per state** — not materialising $m$ rotations of it — so it is
$m/64$ whatever the orbit count turns out to be: **307× at $3^9$, 4,900× at
$3^{13}$**. The orbit count governs *memory* instead, and there $p=7$ is the
case to size for: representatives cost $\text{orbits} \cdot m/8$ bytes against
the current $c\,m^2/8$ for states, so $5^7$'s 2.3 GB becomes kilobytes and even
$p=7$'s $\approx m/2$ orbits stay linear in $m$ where states are quadratic. This
is T1b, and §2.2 is why it is not optional.

**Tests, all shipped and passing (10/10, 51 ms).** (a) Exact agreement with the
tree-condition oracle at every cap, for ten $(p,k,d)$ combinations including
$p=2$ — the fixture the GPU plan's §6 demands. (b) The $k=1$ results reproduce
`drafts/verify_gpu_plan.py` exactly, for both rejected models at once: **cap 0
is the distinct-root model and cap $k$ is the multiplicity model**, so the towers
and states columns at $p = 3,5,7,11,13$ and the landmine counts 78/11/207/371 at
$p = 79,89,97,113$ are asserted verbatim. That is the evidence that the cap
*generalises* both rather than replacing them. (c) The exemplar recovers
$\mathrm{ht}_p$ from its reduced tuple. (d) $V_p(k,d)$ monotone in $k$ and $d$,
towers monotone in the cap. (e) `orbit_misses` is 0. (f) Sizing failures throw
rather than wrap: exact counting past $2^{64}$ and state-table blowup are both
refused.

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
  reported. The floor is $k \ge 11, 7, 5, 4, 3, 3, 2$ at $p = 2, 3, 5, 7, 11,
  13, \le 53$ — and §2.2 now says even *that* is optimistic by a measured factor
  of 2–3, so the working target is $k \approx 15$–$20$ at $p=3$.
- **A depth-4 null at any $k$ is a bug**, since solutions exist — T4 is the
  tripwire, and it has already earned its place: it is what turns "no depth-4
  tower mod $3^5$" from an alarming result into the arithmetic fact $5 < 6$.
- $V_p(k,5)$ growing 1:1 with $k$ **past** the floor, at a prime where
  $V_p(k,4)$ has stabilised, is the signal. Below that it is indistinguishable
  from insufficient precision. **As of the runs in §2 both conditions hold at
  $p = 3$** — $V_3(k,4) = 6$ is stable and $V_3(k,5)$ tracks $k$ at $k = 7, 8$ —
  **and this still is not the signal**, because §2.2's revised floor is not
  cleared. This is exactly the case the rule was written in advance to catch.
- The density branch is a heuristic and stays labelled one, whatever it produces
  (§1's note on the circle method).
- **Report starred entries as starred.** A table that prints $V = k$ without
  marking it invites reading a non-answer as a null; `tower_demo` marks it and
  so should any write-up.
