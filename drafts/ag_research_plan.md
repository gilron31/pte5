# AG research: status and plan (rev. 2)

2026-10-22. This is the working plan for the research branch, separate from the
paper. Evidence is in `ag_research_notes.md` (entries N1–N14), scripts in
`ag_research/`, and the background in `algebraic_geometry_plan.md` (lines G1–G6).

Status labels:

- **PROVEN**: exact computation, or a theorem with proof.
- **EVIDENCE**: strong, but not a proof.
- **NEG**: a negative result within the stated scope.
- **OPEN**: not yet settled.

---

## 1. Where we stand

The two cheapest routes to an $E_5$ through known structure are now closed:

- **Combining Bremner's families:** no two family members form an $E_5$ (PROVEN).
- **Local obstructions:** there are none, at any depth (PROVEN).

The bottom-up route is one step from a first certificate. Symes' $E_4$ has 7
genus-2 quotients of its lift curve with rank exactly 1 (PROVEN), and only the
Sage Chabauty run remains. The one genuinely new arithmetic object is the
exemplar $E_4$ ($L_1 = 33065996717$): it lies on Bremner's slice but on neither
family, and all cheap ways of finding a curve through it have come back empty.

## 2. Steps covered

| line | step | result | status | notes |
| --- | --- | --- | --- | --- |
| G1 | $V_4$ singular locus | 864 ODPs, all degenerate, so $V_4$ is of general type | PROVEN | plan §2.3 |
| G1 | $V_n$ as a Kummer tower; type, degree, $K$ | $K_{V_5} = \mathcal O(18)$ | PROVEN | plan §2.2 |
| G1 | $\mathrm{Sing}(V_5)$ by zero pattern | dimension 1, degenerate: 864 cone lines, 5 curve families, isolated points | PROVEN (stratification) | N12 |
| G1 | local types along $\mathrm{Sing}(V_5)$ | running | OPEN | N12, `v5_local_type.py` |
| G2 | no local obstruction at any depth | $\mathbb Z_p$-towers exist for all $p$ and $d$ | PROVEN | plan §4 |
| G2 | §4 depth-5 floors | settled at $p = 5, 11, 13$; ceilings at 2, 3, 7 | PROVEN (given the §4 lower bounds) | plan §4.3 |
| G2 | top-node local optimality | CMSV is optimal at 5, 13, 17, 37, 47 | PROVEN | N7 |
| G3 | catalog to $10^{12}$ | 15 $E_4$ (3 new beyond $10^{11}$) | PROVEN (certified) | N1 |
| G3 | identities | all starred Table 1 rows are family members; family 1 is closed under twin and is slice $\cap$ hyperplane | PROVEN | N2, N11 |
| G3 | special planes | Bremner's plane is unique; the slice is a $(\mathbb Z/2)^4$-cover of $\mathbb P^2$ over 4 conics with 4 K3 quotients | PROVEN | N6 |
| G3 | the exemplar | on the slice, on neither family, no small relation, isolated through height 1200 | PROVEN / EVIDENCE | N4, N10, N11 |
| G3 | low-degree special curves | none over lines or conics | NEG (proven) | N8b |
| G3 | family-2-type cubic through the exemplar | none, over all 6 lines, coefficients $\le 15$ | NEG (narrow ansatz) | N4a, N4b |
| G4 | $E_5$ from Bremner's families | none: within family 1, within family 2, across | PROVEN (exact) | N8a, N8c, N8d |
| G4 | $\varphi$ on the slice, $\deg 32 = 8 + 24$ | no rational partner with a different $E_4$ in about 34,000 base points or the known points | EVIDENCE | N9 |
| G5 | finiteness per seed | lift curve of genus 98,305 (Faltings) | PROVEN | plan §3.3 |
| G5 | elliptic quotients | all rank $\ge 1$, for structural reasons | PROVEN (Symes) | plan §7.2 |
| G5 | genus-2 descent | implemented; on 56 split Jacobians 0 violations, sharp in 52 | PROVEN (validated) | N13 |
| G5 | rank-1 quotients of Symes' lift curve | 7 certified | PROVEN | N13 |
| G5 | Chabauty + MW sieve on those quotients | all 7 certified: no $E_5$ above Symes' $E_4$, any $\mu$ | PROVEN | N15 |
| G5 | soundness (hide a real point) | the sieve keeps it, 3 of 3 trials | PROVEN | N16 |
| G5 | all catalog $E_4$ + family members (42 seed forms) | 42 of 42 certified: no $E_5$ above any of them, any $\mu$ | PROVEN | N17 |
| G1 | local types of $\mathrm{Sing}(V_5)$ | curves: $cA_1$ (canonical); 128 points: cones over $V_4$ ($a = -3$) | PROVEN (local); $V_5$ general type up to a routine bound | N12a |
| G3 | slice search to height $10^4$ (CUDA) | only the families and the exemplar | EVIDENCE: exemplar sporadic | N18 |
| G4 | $\varphi$ scan, $\lvert t\rvert, u \le 400$ | no slice–slice candidate | EVIDENCE | N18 |

## 3. Action items, by priority

### A. G5: the first bottom-lift certificate (top priority; unblocked by Sage)

- **A1. Smoke-test `chabauty_g2.py`.**
  - Every known rational point must give a zero in its disc.
  - Integrals between two rational points must vanish once $\omega$ is fixed.
  - The tiny integrals must agree with Sage's `tiny_integrals_on_basis`.
- **A2. Run it** on the 7 certified quotients, at the smallest good primes. Where
  a disc's Strassmann bound exceeds its known count, retry at the next prime.
  Several primes can be combined, each disc taking its best bound.
- **A3. Calibrate on a case with a real lift:** the derived $E_3$ of Symes'
  $E_4$, lifting to $E_4$. The fiber has genus 129, and the known lift
  $t = 67405/12$ must survive as a genuine extra rational point. This mirrors
  the synthetic-injection test (§6.8). It needs the descent and Chabauty for 8
  forms, the $E_3$ analog: 4-subsets give elliptic quotients and 6-subsets
  genus 2.
- **A4. Harden the descent.**
  - Add degree-2 divisors (pairs of conjugate local points) to the 2-adic
    sampling, so that no curve is lost to a short image.
  - Fix seeds and log the local images.
  - Add an independent check of a few bounds; Sage's `Jacobian` has none, so
    use Magma if access appears, or a second sampling strategy.
- **A5. Scale up.** For every seed: the 15 catalog $E_4$ in both twin forms, and
  the first $N$ members of each Bremner family (seeds beyond any $B$):
  1. Survey the 6-subsets for rank-bound-1 quotients, about 19% of them for
     Symes.
  2. Certify rank 1 and run Chabauty.
  3. Output per seed: "no $E_5$ above it, for any $\mu$", or "undecided"
     together with the reason.

  Cost: the descent takes about 0.04 s per curve; Chabauty is to be measured.
- **Branch A′ (if a seed gets no certificate):**
  - Try more 6-subsets.
  - Use genus-3 quotients (8-subsets, rank $\le 2$ needed), which needs a
    genus-3 descent.
  - Or run a Mordell–Weil sieve on a rank-1 quotient with a known point.

### B. G1.2: is $V_5$ of general type?

- **B1.** Finish `v5_local_type.py`: is the singular scheme reduced along each
  curve stratum? If so, the transversal singularity is an ODP, which is
  canonical. If primary decomposition stays slow, test at random points over
  $\mathbb F_p$ instead (the local Tjurina number of a transversal slice).
- **B2.** The isolated points: 16 over $K = (0)$, and the 128 points where all 8
  branch divisors meet. Take local equations and compute log discrepancies by
  one blow-up (or the Newton polyhedron), and decide canonical or not.
- **B3.** If everything is canonical, then **$V_5$ is of general type**
  (Theorem). If not, estimate the plurigenera; $K = \mathcal O(18)$ leaves room.

### C. Q1: is the exemplar on a third special curve?

- **C1. Push the slice search to height about 5000.** Use a C++ port of the
  QR-sieve, or a sieve in the $(s:t:u)$-model. A third rank-1 curve with
  moderate height growth should produce companions there.
  - **Gate:** 2 or more new off-family slice points: fit curves through them,
    which could be a third family. Nothing up to 5000: treat the exemplar as a
    sporadic point.
- **C2. Use the K3 quotients.** For $K_{bdf}$ and the other three, take the
  elliptic fibrations (pencils of lines through the tacnodes) through the
  exemplar's points, compute the fiber ranks in PARI, and intersect the
  candidate curves on the four K3s.
- **C3. Picard number of $X_1$** (needs a C++ point counter over
  $\mathbb F_{p^n}$). That gives the lattice of split curves and a systematic
  list of low-degree candidates.
- **C4. A broader ansatz:**
  - rational quadrics $Q$, allowing denominators;
  - family-1 type, with $T = F + H$ a cubic form of a different shape;
  - quartic curves.

### D. $E_5$ with both halves on the slice ($\varphi$-correspondence)

- **D1. Geometry of the degree-24 correspondence $X'$.** Determine its genus over
  a line, its irreducibility (the generic fiber already is), and its Kodaira
  dimension. That decides whether rational points of $X'$ can come in families.
- **D2. A larger scan:** $\lvert t\rvert, u \le 1000$, plus every slice point
  found in C1.

### E. Local theory (G2) follow-ups

- **E1. Better witnesses from $\mathbb Z_p$-seeds,** not only integer ones (the
  same trick, applied one level down). The aim is lower ceilings at
  $p = 2, 3, 7, 11$ and the mid-size primes.
- **E2. Regression tests for the C++ local model:** $V_3(15,5) \le 14$,
  $V_7(11,5) \le 10$, $V_{11}(5,5) = 4$. They become testable once the orbit
  rewrite reaches those precisions.

### F. Hygiene

- **F1.** Write a `ag_research/README.md`: a script index, the run order,
  expected outputs and runtimes.
- **F2.** Keep the operational rules:
  - no `modStd`, which forks a worker per core regardless of `setcores`;
  - timeouts on Singular jobs;
  - no `pkill -f` with a pattern contained in its own command line;
  - in PARI, set `parisizemax` in the top-level script only, never in a file
    passed to `read()`.

## 4. Order of work

1. **Now:** A1–A2 as soon as Sage imports; B1 runs in the background.
2. **Next:** A3 (calibration), A4, then A5 on all catalog seeds.
3. **In parallel, cheap:** C1 (C++ sieve), D2, E1.
4. **Later:** B2–B3, C2–C3, D1.

## 5. What would change the plan

| outcome | consequence |
| --- | --- |
| A2 certifies Symes | G5 becomes a pipeline (A5); its output is a list of seeds with no $E_5$ above them for any $\mu$ |
| A2 fails at every prime tried | there are rational points beyond the known 8: an $E_5$ above Symes, or a bug. Check A1/A3 first |
| C1 finds new off-family slice points | Q1 becomes "find the third family", which is the G3 jackpot; feed its members into G4 and A5 |
| B2 finds non-canonical points | general type for $V_5$ needs a plurigenus estimate; Lang–Vojta claims about $V_5$ must wait |
| D1 shows $X'$ is not of general type | slice–slice $E_5$ candidates could come in families; scan along the special curves of $X'$ |
