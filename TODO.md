# TODO

Target: the fused "foundations + experiments" paper, organized **bounds-forward**
(E₅ lower bounds as the headline), committing to the full exhaustive campaign.
Section skeleton and inline `**TODO:**` gaps live in `paper.md`; the structural
brainstorm and rationale are in `~/.claude/plans/read-the-entire-repo-zany-kitten.md`.
Detailed drafts of most items below are in `drafts/critical_review.md`.

Sequencing: **W0 (literature) first** — it changes what W1/W2 may claim — then
W1 (theory) and W2 (code) in parallel, then the campaign, then W3 write-up.

## W0 — Literature pass (BLOCKING, do first)

- [ ] Identify reference [3] from `sources/Bre08.pdf`; extract both families in
      symbolic form $L_1(k),L_2(k),L_3(k)$ (blocks the family-intersection run).
- [ ] Read `sources/CSMV23.pdf` (CMSV) in full, esp. §3 local families —
      decides whether "sharper required divisors" (§9) is a theorem or foreclosed.
- [ ] Read `sources/BMR09.pdf` in full — do the Litter conditions already
      articulate the tower recursion? Governs every "we introduce" in §3.
- [ ] Classical sweep: Chernick 1937, Gloden 1944, Dickson vol. 2,
      Dorwart–Brown, Caley's ℤ[i] thesis, Choudhry, Borwein–Lisoněk–Percival.
- [ ] Existing tables (Chen Shuwen, Meyrignac eslp) — overlap with the catalog?
- [ ] Cite `sources/TLT26.pdf` (Tsai et al., E₃-density) in §8.

## W1 — Theory (paper §§2–5, 8)

- [ ] Easy: §2 corrected equivalence + Tarry counterexample; "recursively
      symmetric subclass" framing.
- [ ] Medium: §3 tower recursion theorem (uₖ = (r²₂ₖ − r²₂ₖ₊₁)/2) with explicit
      2-power bookkeeping + parity corollary + bottom-lift criterion.
- [x] §4 divisibility ladder (Thm 4.1) + E₅ lower-bound corollary (4.2) — DRAFTED.
      Constants independently re-derived and **corrected**: √L₅ ≥ C′₁₆/2 ≈ 2.44×10²⁷
      (the review's 1.43×10²⁶ dropped the prime 17); root > 850, L₁ > 3.6×10⁵ from
      CMSV's height-850 search (the ladder alone gives only L₁ ≳ 2650). Remaining:
      confirm C′₁₆ exponents directly vs CMSV Table 1.
- [x] §5 **surjectivity for all admissible norms** (repeated primes, units,
      ramified 2) — the load-bearing lemma. DONE: Prop 5.1 (constructive token
      argument), verified `verify_claims.py` §6.
- [x] §5 **single-constraint reduction** $E_4 \iff W=0$. DONE: Prop 5.2; with
      Prop 5.1 gives "$E_4$ solutions = ℤ[i]-points of $\{W=0\}$".
- [x] §5 degenerate locus of W — DONE (Prop 5.4: genuine E₄ ⟺ W=0 ∧ Δ≠0, where
      Δ = the coordinate discriminant; D = {W=0}∩{Δ=0} is proper closed).
- [ ] Hard: §5 canonization lemma (unique syndrome representative; the b≥2
      multiset case) — still open, blocks the §7 census.
- [x] §5 ℤ[√2] theorem (L₂,d odd primes ≡ ±1 mod 8) + smoothness — DRAFTED (Thm 5.3).
- [ ] §8 density heuristic — research PLAN drafted (multidegree H⁻²/H⁻¹⁸, r₂
      second-moment lift of [TLT], 5-step plan). Computation pending: confirm V₅
      multidegree/canonical class; σ_p singular series (needs §6.6 Fₚ run); the r₂
      correlation-sum main term for E₄.
- [ ] Misc: trace the etymology of the "Litter" conditions.

## W2 — Code / campaign (paper §§6–7)

- [ ] Verification harness: tower recovery, degree-16 polynomial-identity check,
      ladder-divisibility check for every cataloged solution (correctness spine).
- [ ] Extend the enumerator to **non-squarefree norms** (the S3 hole); every
      exclusion claim rests on this. Math now discharged (Prop 5.1); implement the
      exponent-vector iteration per its proof — reference recipe in
      `verify_claims.py` §6 (`reconstruct`).
- [ ] Family self-test: re-find all of [3]'s members incl. non-squarefree k.
- [ ] Meet-in-the-middle in the extended form, in C/C++/CUDA.
- [ ] Exhaustive campaign, two pre-registered regimes: (i) exhaustive over
      canonical L₁ ≤ B (B from pilot, 10¹⁰–10¹²); (ii) champion norms (~20 split
      primes) on GPU. Keep separated in the write-up.
- [ ] Dedup mod the full symmetry group + syndrome census via canonization.
- [ ] Top-collision search on (L₁,L₂,L₃); bottom-lift screen (μ to 10⁵–10⁶,
      parity wheel); Fₚ full-splitting tower sieve (p ≤ 10³, depths 3–5).
- [ ] Symbolic family intersection (resultants/Gröbner) — blocked on W0 [3].
- [ ] Hardening protocol: two independent implementations, synthetic injection,
      checkpointed logs, Fₚ count cross-checks.
- [ ] Statistics package for §7.

## W3 — Writing / infra

- [x] §6 rewritten as a full search design (6.1–6.8, incl. Prop 6.1 exclusion);
      §8 drafted as a density research plan. (Prose only — the runs are W2.)
- [x] Defect fixes (critical pass, see `drafts/sources_insight.md` §5c): §2
      level-3 tree constant **32L₃ → 128L₃** (verified; general constant 2^{2^k-1});
      abstract correspondence (single E_n → size 2^{n-1}, not "in pairs"); §1
      τ-conjecture implication direction; §4 required-divisor definition; trimmed
      §4 internal "correction" note.
- [ ] Rewrite the abstract bounds-forward.
- [ ] Write §1 Introduction (PTE state of the art, CMSV, BMR/gems, TLT, contributions).
- [ ] Fix references: full citations, identify [3], add [4] CMSV, [5] TLT.
- [ ] Fill the Appendix exemplar tower (verify via harness before citing).
- [ ] Settle scope/credit among authors before heavy drafting (ladder,
      recursion, bounds originated in the review, not the original proposal).
