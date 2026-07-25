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
- [ ] Medium: §4 divisibility ladder theorem + E₅ lower-bound corollary;
      **independently re-derive** the constants (√L₅ ≥ 1.43×10²⁶, root > 3400,
      L₁ > 5.8×10⁶) from CMSV's C′₁₆.
- [ ] Hard: §5 **surjectivity for all admissible norms** (repeated primes,
      units, ramified 2) — the load-bearing lemma; + canonization lemma; +
      degenerate locus of W.
- [ ] Medium: §5 ℤ[√2] theorem (L₂,d odd primes ≡ ±1 mod 8) + smoothness
      explanation.
- [ ] Hard: §8 density heuristic / multidegree computation; geometry framing.
- [ ] Misc: trace the etymology of the "Litter" conditions.

## W2 — Code / campaign (paper §§6–7)

- [ ] Verification harness: tower recovery, degree-16 polynomial-identity check,
      ladder-divisibility check for every cataloged solution (correctness spine).
- [ ] Extend the enumerator to **non-squarefree norms** (the S3 hole); every
      exclusion claim rests on this.
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

- [ ] Rewrite the abstract bounds-forward.
- [ ] Write §1 Introduction (PTE state of the art, CMSV, BMR/gems, TLT, contributions).
- [ ] Fix references: full citations, identify [3], add [4] CMSV, [5] TLT.
- [ ] Fill the Appendix exemplar tower (verify via harness before citing).
- [ ] Settle scope/credit among authors before heavy drafting (ladder,
      recursion, bounds originated in the review, not the original proposal).
