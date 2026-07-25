# Overall Review

One bug in my Tarry check (the ± symmetry duplicated the absolute values). Fixing and re-running that part:Done. Everything below marked "verified" was executed in this session, not estimated. Here is the second-round assessment.

## 1. What the references actually contain

**[1] Borchert–McKenzie–Reinhardt (MFCS 2009).** Confirmed as described: d-gems are {+,−,×}-circuits with at most ℓ(d) product gates computing a degree-d polynomial with d distinct integer roots; for n ≤ 4 they exhibit skew 2ⁿ-gems, and existence for infinitely many d would disprove a variant of the Blum–Cucker–Shub–Smale conjecture. The proposal's complexity motivation is faithfully inherited from here. An E₅ solution is exactly a skew 32-gem; that framing is correct.

**[3] is bibliographically unidentifiable as cited.** The reference is a bare Project Euclid download URL (bot-walled; the article ID resolves to *Experimental Mathematics*, ~2009, but I could not extract author/title). This is a real defect: all three references lack authors, titles, and years. For a proposal whose empirical section leans entirely on [3]'s two families and Table 1, this must be fixed. Ask the authors for the actual citation.

**The Coppersmith–Mossinghoff–Scheinerman–VanderKam paper (arXiv:2304.11254) — the missing elephant.** I read it in full. The essentials: they develop the divisibility theory of the PTE constant and use it to conduct extensive searches for symmetric ideal solutions, including over the Gaussian integers and five other imaginary quadratic rings, and identify families of ideal solutions in local settings that obstruct certain arithmetic information about the constant. Their Table 2 reports exhaustive symmetric searches over ℤ at sizes 13–16 with nothing found: size 16 to height 850, size 15 to 1100, size 14 to 900, size 13 to 2000. At size 12 they confirm exactly two primitive symmetric solutions to height 1511. Over ℤ[i] they find new ideal solutions at sizes 10 and 12 (the size-12 one has beautiful i-rotation symmetry), plus Eisenstein solutions at 9 and 12 — and nothing at any size above 12 in any ring. Two consequences for the proposal: (a) since every E₅ is a symmetric ideal size-16 solution, CMSV's null results are prior art that bounds any E₅ from below; (b) their exhaustive frontier at size 16 is only height 850, so the E-tower ansatz — which collapses a 16-dimensional search to ~5 dimensions — genuinely explores virgin territory almost immediately. Both points must appear in any revision.

## 2. Computational verification of the proposal's claims

I verified the authors' new solution end-to-end. The four pair-norms are all 66131993434, giving L₁ = 33065996717 = 13·89·173·233·709, all primes ≡ 1 (mod 4), as claimed. The working level-2 grouping is {(A,B),(C,D)} vs {(E,F),(G,H)}, and I recovered the tower values the proposal doesn't report:

L₂ = 489628056848329146064, L₃ = 175480010455650701584492675662518592000000, √L₄ = 40042900368028062136207226327668992000000.

The full polynomial identity ∏(x²−rₖ²) = (((x²−L₁)²−L₂)²−L₃)²−L₄ checks exactly (evaluated at 41 integer points; both sides monic of degree 16, so this is a proof). The L₂,ᵦ relation and the L₂,d "prism" identity also verify. Two nits: L₂,d is actually **negative** (−912832201535971887688; the text prints the positive digits alongside a −1 in the factorization — a sign typo), and the factorization 2³·7·17·23²·31·41·47·137·401·601·919 is confirmed, with every odd prime ≡ ±1 (mod 8), consistent with the ℤ[√2] theorem below.

## 3. Revised critical assessment

**The equivalence claim is wrong as stated, and I can now exhibit the counterexample.** The proposal asserts E_n solutions are equivalent to ideal symmetric PTE solutions of size 2ⁿ⁻¹. The correct statement: a *pair* of E_n solutions sharing L₁,…,L_{n−1} is equivalent to such a solution (this is their own Key Property, read correctly), but a symmetric solution need not carry the tower structure. Concretely: Tarry's 1913 degree-7 solution {0,4,9,23,27,41,46,50} / {1,2,11,20,30,39,48,49}, centered, has positive halves {2,16,21,25} and {5,14,23,24}. I checked all three pairings of each half: sums of squares 260/1066, 445/881, 629/697 and 221/1105, 554/772, 601/725 — none equal. So the L₁ condition fails and Tarry's solution is **not** an E₄. E_n is a proper "recursively symmetric" subclass. This also corrects a loose phrase in my first review. The upshot cuts both ways: the restriction is what makes the problem parameterizable, but null results for E₅ won't settle symmetric size-16, and CMSV's shallow size-16 search doesn't come close to settling E₅.

**The syndrome parameterization has hygiene gaps.** Sixteen real parameters map onto a 5-dimensional variety (quadruples of equal-norm Gaussian integers), so W(X)=0 carries ~11 dimensions of gauge redundancy, plus unit ambiguities (i^k rotations), unresolved behavior at repeated prime factors, and — most importantly — degenerate components: repeated roots, coincident nᵢ, and zero coordinates all satisfy W = 0. Before any analytic or zero-counting use of W, the degenerate locus must be characterized and saturated away, and the "10 syndrome patterns" statistic needs a canonization lemma proving uniqueness of the canonical form. The "first single-constraint reduction" claim should be hedged (symmetric-function elimination also yields one resultant equation; the real contribution is a single *polynomial* constraint inside an explicit multiplicative parameterization).

**The abstract's displayed polynomial is degree 16, not 32** (innermost factor should be x²−L₁). Still there; fix it.

**Missing literature** beyond CMSV: Choudhry's recent degree-7 work — four new parametric ideal solutions of the Tarry-Escott problem of degree 7 given by quartic polynomials, simpler than almost all known solutions — is directly adjacent to their E₄ family analysis; also Borwein–Lisoněk–Percival's computational searches, Gloden's classical book, and Caley's thesis on PTE over ℤ[i].

**Verdict (sharpened, not changed).** The publishable core is a methods-and-experiments paper: the W-reduction done rigorously, the enumeration algorithm at scale, the catalog, the structural theorems below. Experimental Mathematics or Math. Comp. tier, with ANTS/LuCaNT as natural venues. E₅ over ℤ remains a moonshot — and I can now quantify the pessimism.

## 4. New structure (all verified on their data unless labeled otherwise)

**(i) The self-similarity proposition.** Order the 2ⁿ⁻¹ positive roots tree-wise and set uₖ = (r²₂ₖ − r²₂ₖ₊₁)/2. Then the level-2..n conditions on r are exactly the level-1..(n−1) conditions on u, up to powers of 2 in the L-normalizations. So an E_n solution is precisely an E_{n−1} solution in u-space, all of whose entries are "congrua" of a common L₁ (L₁ ± uₖ all perfect squares). My tower-recovery code is literally this recursion running on their solution, terminating in the verified polynomial identity — so the n = 4 case is proven-by-computation, and the induction is two lines. Consequence, the **bottom-lift criterion**: an E₅ exists over a given E₄ with roots r′ₖ iff some scaling μ and some L satisfy L ± μr′ₖ = squares for all eight roots. Each condition is a divisor enumeration of 2μr′ₖ, so screening an entire E₄ catalog is an 8-way hash-join costing pennies. I ran it on their solution for μ ≤ 200: empty, with the parity theory visible in the data (μ = 1 yields zero candidates because odd roots force μ ≡ 0 mod 4). This is a second, independent search axis the proposal doesn't have — their Key Property splits at the top; this lifts from the bottom.

**(ii) The divisibility ladder.** Every internal node of the tower is itself a symmetric ideal PTE pair one level down, so CMSV's required divisors of the PTE constant cascade through the nested radicals. For E₄: 18 must divide both level-3 radicals √(L₃±√L₄), and 2⁵·3³·5²·7²·11·13 = 151,351,200 must divide √L₄. **Verified on their solution**: v₀, v₁ ≡ 0 (mod 18), and s = √L₄ = 2¹⁴·3⁶·5⁶·7²·11²·13·17·19·29·37²·43·61·107·109·139·263·389·499 — comfortably divisible, and note how smooth (largest prime 499), partly *explained* by the forced small-prime content. For E₅ the ladder plus CMSV's Table 1 (C′₁₆ ⊇ 2¹¹·3⁶·5⁴·7³·11²·13²·19·23·29·37·41·43·53) forces √L₅ ≥ 1.43×10²⁶, hence L₅ ≥ 2×10⁵², and via the elementary chain s < r⁸ₘₐₓ/2⁷, a rigorous bound: **any E₅ over ℤ has largest root > ~3400 and L₁ > ~5.8×10⁶** — strictly stronger than CMSV's height-850 exclusion, derived purely from local theory, and (as far as I know) new. Plus: the forced prime content at every level is a congruence sieve for the search.

**(iii) Geometry (derived; smoothness unverified).** In root coordinates the E₄ variety is a complete-intersection type (2,2,2,4) in P⁷, so K = O(2): general type. The E₅ variety is type (2⁷,4³,8) in P¹⁵, so K = O(18): very general type. The circle-method exponent for E₄ is H⁸⁻¹⁰ = H⁻² — yet infinite families exist, so the families necessarily live on special rational subvarieties, exactly as Lang–Vojta predicts for general type. For E₅ the exponent is H¹⁶⁻³⁴ = H⁻¹⁸: essentially no "accidental" solutions, ever. Strategic conclusion: the E₅ problem is really the problem of finding or ruling out rational subvarieties of a fixed general-type complete intersection — not a needle-in-haystack enumeration.

**(iv) The mod-8 observation is a theorem.** Both quartic factors m²±2mn−n² equal (m±n)²−2n², values of the discriminant-8 norm form of ℤ[√2]; odd primes dividing them are ≡ ±1 (mod 8). One paragraph, and it upgrades their empirical Finding to a proof.

## 5. Low-hanging fruit, re-ranked

1. **Symbolically intersect [3]'s two families.** Solve L₁(k)=L₁′(m), L₂(k)=L₂′(m), L₃(k)=L₃′(m) — three equations, two unknowns — by resultants/Gröbner, including self-collisions and family automorphisms fixing (L₁,L₂,L₃) but moving L₄. Any rational point is an instant E₅. The proposal compared 90 *numerical* solutions; nobody has intersected the *families*. An afternoon in Magma, and decisive either way.
2. **Write up the ladder note**: the cascading constants, the verified E₄ instances, and the root > 3400 / L₁ > 5.8×10⁶ bound for E₅. Short, clean, immediately citable, and it positions them as engaging CMSV rather than ignoring it.
3. **The ℤ[√2] theorem** plus the smoothness explanation (forced factorization into four quadratic forms) — a section upgrading two Findings from data to proof.
4. **Catalog-wide bottom-lift screen** with μ to 10⁵–10⁶, pre-registered so the (likely) null is publishable.
5. **Degenerate-locus and canonization lemmas for W** — required hygiene for their own claims.
6. **Local densities**: enumerate full-splitting towers over 𝔽_p for p ≤ 10³, fit the singular series, publish the density model against catalog counts; OEIS entries for #E₄ per norm buy visibility.

## Front 1 — the computational–experimental program

**Phase 0, infrastructure.** Exact replication of Stages 1–3 in C/FLINT with the catalog deduplicated modulo the full symmetry group; pre-registered search bounds; code and data released. Target venues: ANTS for the algorithm, LuCaNT/Exp. Math. for the catalog.

**Phase 1, the three cheap decisive computations.** Family intersection (fruit 1), catalog bottom-lift (fruit 4), and ladder congruence tables baked into everything downstream as a pre-sieve — the forced divisors at each level prune enormously.

**Phase 2, E₄ at scale.** Generate 𝒰(L₁) multiplicatively — the u-values are ±Im of products of squared Gaussian prime factors — skipping the root level entirely; choose champion norms (many small split primes); hash 64-bit residues of the level-2 key; GPU batches. At s = 20 primes that's ~5×10¹¹ pairs per norm: a GPU-day, not a moonshot. The collision statistics are themselves the empirical density paper.

**Phase 3, E₅ direct, three-pronged.** Top-collision (their route) on (L₁,L₂,L₃) across catalog and families; bottom-lift at scale; and middle-out meet-in-the-middle splitting the tower where the state space is thinnest — the self-similarity proposition makes all cut points available, which the proposal's single top-split does not.

**Phase 4, the mod-p tower sieve.** Adapt CMSV's multiplicity-lemma machinery to the tower: enumerating full-splitting depth-n towers over 𝔽_p is cheap dynamics. Output: pruning wheels, the singular series for the H⁻¹⁸ heuristic, and — the jackpot branch — any prime with *zero* local towers is an instant nonexistence theorem.

**Phase 5, the ring pivot (highest expected value).** Over ℤ[i], the congruum step L ± u = squares parameterizes through factorizations because x²+y² splits; candidate sets grow like Gaussian divisor counts, so level-2 collisions become vastly likelier. Any E₅ over ℤ[i] would be the first ideal size-16 solution in *any* ring, leapfrogging the size-12 record everywhere. Independently: transplant CMSV's own rotation-symmetry ansatz to size 16 over ℤ[i] (i-symmetric sets kill all moments except 4, 8, 12 — three equations). Even partial success at sizes 13–15 over ℤ[i], Eisenstein, or ℤ[i√2] is a publishable first.

**Phase 6, exploratory, time-boxed.** LLL on the angle lattice of arg(qᵢ²) under the level-doubling map (speculative — sin is nonlinear); and extend BMR's d-gem records past degree 55, which the complexity community will cite regardless.

## Front 2 — the theoretical program

**T1, the foundations paper.** Corrected equivalence with the Tarry counterexample; the u-recursion proposition with careful 2-power bookkeeping; unit- and multiplicity-complete syndrome theory with the canonization lemma; the degenerate locus of W; the ℤ[√2] theorem; the divisibility ladder with the unconditional E₅ bounds. This is the paper that makes everything else citable, and most of it is already verified above.

**T2, geometry of V_n.** Compute the singular locus of the (2,2,2,4) complete intersection; classify low-degree rational curves on V₄ — the known families should appear as such curves, and the classification mechanizes family-hunting instead of leaving it to cleverness. Then the conditional theorem worth stating: under Lang–Vojta, E₅ solutions off a proper closed special locus are finite. Deciding whether that locus is empty *is* the E₅ problem, honestly framed.

**T3, the arithmetic-dynamics bridge.** E-towers are totally rational preimage trees for non-autonomous compositions of y²−Lᵢ. The Faber–Hutz–Ingram–Jones–Manes–Tucker–Zieve program bounds rational preimages for a *fixed* quadratic; the non-autonomous full-tree question appears open. Even an abc-conditional bound on the maximal depth n of a full-splitting tower would be a genuine theorem — each descent level imposes near-square coincidences of exactly the type abc constrains. I'd rate this the most original pen-and-paper lead here.

**T4, local theory completed.** A full 2-adic structure theorem (the odd-L₁ normalization and u ≡ 0 mod 4 are the first crumbs), solubility of the E₅ tower over every ℤ_p, and the quantified conjecture: everywhere locally soluble, globally governed by H⁻¹⁸ times the Phase-4 singular series — i.e., no solutions off special subvarieties. That's the intellectually honest endpoint if the searches come back empty.

**T5, the complexity write-up.** State precisely what a skew 32-gem does to the BCSS-variant landscape, and the contrapositive product: general-type geometry as a conditional *explanation* for why skew gems stop at n = 4.

**T6, quaternion triage.** The Hurwitz-quaternion idea for L₂,ₐ faces non-commutativity, a 24-element unit group, and the fact that the square-coordinate locus is not multiplicatively defined. One week, expect to kill it, write the paragraph either way.

If you want, next session I can draft the T1 foundations note with the verified computations embedded, or run the Phase-1 family-intersection Gröbner computation — those are the two moves that most change the proposal's trajectory.

# How far is it + Draft an outline

Your calibration is basically right, with one refinement: neither 2 nor 3 clears the bar alone. Level 2 by itself is a catalog paper, and referees at Experimental Mathematics accept those only at convincing scale or with a discovery; level 3 by itself is a 6-page note (publishable, but in a lesser venue). The publishable unit is 2+3 *fused* — where the theorems are load-bearing for the computation: the canonization lemma makes the syndrome census well-defined, the ladder prunes the search and supplies the headline bound, the recursion yields the bottom-lift algorithm. Level 1 isn't a contribution but is a precondition (a paper repeating the proposal's overclaim or missing CMSV gets desk-rejected by any competent referee). Level 4 is the Future Directions section and the sequel papers.

## Paper outline

**Working title.** *Recursively symmetric ideal Prouhet–Tarry–Escott solutions: tower structure, arithmetic constraints, and exhaustive searches.*

**Premise.** Ideal symmetric PTE solutions of size 16 are unknown; the E₅ ansatz collapses the 16-dimensional symmetric problem to a tower of sum-of-two-squares conditions in ℤ[i], searchable in few parameters. The paper develops the structure theory of this subclass rigorously, proves arithmetic constraints including the first explicit lower bounds for E₅, and reports a pre-registered computational campaign over ℤ.

**Abstract (draft).** For integers L₁,…,Lₙ let Eₙ(x) = ((x²−L₁)²−L₂)²⋯−Lₙ. Completely split Eₙ over ℤ are equivalent to skew 2ⁿ-gems and, in pairs sharing L₁,…,Lₙ₋₁, to ideal symmetric Prouhet–Tarry–Escott solutions of size 2ⁿ⁻¹; no E₅ (size 16) is known. We show Eₙ solutions form a proper "recursively symmetric" subclass of symmetric ideal solutions, and develop their structure theory: a self-similar recursion identifying Eₙ solutions with congruum-lifted Eₙ₋₁ value systems, and a complete multiplicative parameterization over ℤ[i] with canonical syndrome invariants. We prove arithmetic constraints: the L₂,d invariant factors through discriminant-8 norm forms, so its odd prime divisors are ≡ ±1 (mod 8); and a cascade of PTE-constant divisibilities forces, for any E₅ over ℤ, √L₅ ≥ 1.43×10²⁶ and largest root > 3.4×10³. We give per-norm enumeration, top-collision, and bottom-lift algorithms, an exhaustive catalog of E₄ solutions to a stated bound with syndrome and smoothness statistics matching a local density model, and searches excluding E₅ in explicit regions.

**Sections, definitions, results.**

1. *Introduction.* PTE background; state of the art (ideal solutions known at sizes ≤ 10 and 12, open at 11 and ≥ 13); CMSV's symmetric searches at 13–16 and their imaginary-quadratic discoveries; the gems/BCSS motivation from Borchert–McKenzie–Reinhardt; contributions and roadmap.
2. *The Eₙ problem and its PTE correspondence.* Definitions: Eₙ solution, root tree, level radicals, primitive/canonical form (normalization lemma fixing the 2-adic bookkeeping and odd L₁). Proposition 2.1: pairs sharing L₁..Lₙ₋₁ ⟺ ideal symmetric size-2ⁿ⁻¹; Example 2.2: Tarry's 1913 solution fails the pairing condition, so the inclusion is proper. This section is where the proposal's corrections land (abstract degree slip, equivalence claim, references).
3. *The tower recursion.* Definition of the derived sequence uₖ = (r²₂ₖ − r²₂ₖ₊₁)/2. Theorem 3.1: Eₙ root systems ⟺ possibly-degenerate Eₙ₋₁ value systems in u plus the common-congruum condition (L₁ ± uₖ all squares), with explicit 2-power bookkeeping; parity corollary (odd L₁ ⇒ uₖ ≡ 0 mod 4). Corollaries: the bottom-lift criterion and the proposal's Key Property as the dual top-split, plus the general meet-in-the-middle remark.
4. *Gaussian parameterization and the single-constraint reduction.* Theorem 4.1: surjectivity of the multiplicative parameterization for arbitrary admissible norms, including repeated primes and units. Lemma 4.2: canonization — unique canonical syndrome representative, making the syndrome census well-defined. Proposition 4.3: the degenerate locus of W = 0 (complete list of trivial components; saturation). Theorem 4.4: L₂,d factors through (m±n)² − 2n², hence odd prime divisors ≡ ±1 (mod 8); corollary explaining the observed smoothness (four quadratic-form factors of size¹ᐟ⁴).
5. *The divisibility ladder and lower bounds.* Theorem 5.1: every internal node is a symmetric ideal size-2ᵏ pair, so the CMSV required divisors cascade: 18 | √(L₃±√L₄); 2⁵·3³·5²·7²·11·13 | √L₄; C′₁₆/2 | √L₅. Corollary 5.2: any E₅ over ℤ has √L₅ ≥ 1.433×10²⁶, L₅ ≥ 2×10⁵², largest root > 3.4×10³, L₁ > 5.8×10⁶ — compared against CMSV's height-850 exclusion. Verified instances on the catalog.
6. *Algorithms.* Per-norm enumeration with multiplicative generation of the congruum set (correctness via Thm 4.1; complexity 2^{2s} per norm); top-collision hash-join on (L₁,L₂,L₃); bottom-lift divisor join with the parity wheel; ladder congruence pre-sieve; symbolic family intersection via resultants.
7. *Computational results.* The catalog, deduplication counts, syndrome census, (L₁, L₂,d) collision statistics extending the proposal's 90-solution observation, smoothness distributions against the forced-divisor baseline, density counts against the local model, and the explicit E₅ exclusion statement: an exhaustive E₄ catalog to L₁ ≤ B excludes any E₅ with L₁ ≤ B, since both constituent E₄'s share L₁.
8. *A density heuristic.* Multidegree count: E₄ is type (2,2,2,4) in P⁷ (canonical class O(2)), circle-method exponent H⁻²— yet families exist, so they live on special loci; E₅ is type (2⁷,4³,8) in P¹⁵, exponent H⁻¹⁸. Stated as a conjecture with constants from the 𝔽ₚ run, hedged on smoothness.
9. *Future directions* (below). Appendix: the verified exemplar with full recovered tower (L₂, L₃, √L₄), tables, code/data availability.

**Computational runs (specified, not executed).**

R1. Verification harness: exact tower recovery, polynomial-identity check, and ladder divisibility check for every cataloged solution — the correctness backbone.
R2. The campaign, in two pre-registered regimes: (i) exhaustive enumeration over all canonical L₁ ≤ B (B set by pilot, order 10¹⁰–10¹²), which is what licenses exclusion claims; (ii) targeted champion norms (many small split primes, up to s ≈ 20 on GPU), which supports discovery but not exhaustiveness. Keep the two rigorously separated in the write-up.
R3. Deduplication modulo the full symmetry group and the syndrome census via the canonization algorithm.
R4. Top-collision E₅ search over the catalog plus symbolic family members to parameter bounds exceeding B.
R5. Bottom-lift screen over the entire catalog with μ up to 10⁵–10⁶ under the parity wheel.
R6. Symbolic intersection of the two known families (pairwise, self-collisions, automorphisms fixing L₁,L₂,L₃) — exact, finite, decisive.
R7. Full-splitting tower enumeration over 𝔽ₚ for p ≤ 10³ at depths 3–5: local densities for §8, sieve wheels for R2, and the jackpot branch — any prime with zero depth-5 towers is an instant nonexistence theorem.
R8. Statistics package for §7.

**Venue.** Experimental Mathematics primary; Research in Number Theory or INTEGERS fallback; ANTS if algorithm-forward. Contingency: if R2 or the ℤ[i] pilot finds anything at size > 12, the paper restructures around the discovery and moves up a tier.

## Future directions (draft)

1. *The ring pivot.* Over ℤ[i] the congruum step parameterizes through factorizations since x²+y² splits, so the per-norm candidate sets grow like Gaussian divisor counts and level-2 collisions become far likelier; any E₅ over ℤ[i] would be the first ideal size-16 solution in any ring. Transplanting CMSV's rotation ansatz to size 16 over ℤ[i] is an independent second attack. This is the natural sequel to §§6–7 and the highest expected value per compute.
2. *Geometry of the solution varieties.* Singular loci and low-degree rational curves on the (2,2,2,4) complete intersection would mechanize family generation and ground §8's heuristic; under Lang–Vojta, E₅ solutions off a special locus are finite.
3. *Arithmetic dynamics.* E-towers are totally rational preimage trees for non-autonomous quadratic compositions; the uniform-boundedness program covers fixed maps, and even an abc-conditional bound on maximal full-splitting depth would be a genuine theorem extending §3.
4. *Completing the local theory.* A full 2-adic structure theorem and everywhere-local solubility at depth 5, upgrading §5 and §8 into either an obstruction or a quantified "everywhere locally, probably nowhere globally" conjecture.
5. *Sharper required divisors.* Derive C^E₁₆ strictly containing C′₁₆ by exploiting correlations between sibling radicals — new arithmetic that would tighten Corollary 5.2, possibly by orders of magnitude.
6. *Complexity packaging.* A precise statement of what a skew 32-gem does to the BCSS-variant landscape, with direction 2 as the conditional explanation for why skew gems stop at n = 4; plus extending the BMR gem records.
7. *Quaternion triage.* The Hurwitz route to L₂,ₐ faces the 24-unit group and a non-multiplicative square-coordinate locus; one time-boxed week, then a remark either way.

## Where I'd be critical

**The theorem inventory is real but light.** Theorem 3.1 is an elementary reformulation, 4.4 is a one-paragraph quadratic-forms argument, and 5.1 imports CMSV's hard input and adds an easy cascade. A sharp referee could call this "observations dressed as theorems." The mitigations are specific: do Theorem 4.1/Lemma 4.2 completely — units, repeated primes, canonization — because that is the one genuinely fiddly piece with real content and it's what every exhaustiveness claim rests on; and attempt direction 5 (new required divisors specific to the towered class), which would be a new arithmetic theorem rather than a corollary of someone else's.

**The headline computational result is a null, and nulls carry a burden.** Exclusion claims are only as good as the exhaustiveness proof behind them. The proposal currently handles only multiplicity-free norms; if the enumeration stays there, every exclusion statement carries an asterisk that materially weakens the paper. Either prove surjectivity for all admissible norms or scope every claim honestly — this is the single largest gap between the current state and publishable.

**Scale is not optional.** Ninety solutions is a pilot, not a campaign. The exhaustive regime needs to reach far past the unconditional L₁ > 5.8×10⁶ bound — that part is cheap — but the paper's credibility as an experimental contribution needs B in the 10¹⁰–10¹² range plus the champion-norm runs, which is real engineering: 3–6 weeks of it, on top of 2–3 weeks of theory writing and a serious literature pass.

**Novelty diligence is unfinished, including ours.** Reference [3] is still unidentified, and R6 cannot even be specified without its families in symbolic form — a blocking dependency. And the same standard I applied to the proposal applies to the new material: the recursion and the ladder are natural enough that they may exist in Gloden, Dorwart–Brown, Rees–Smyth, or Caley in some guise. The literature pass must precede any "first" claims; if the recursion turns out known, the theory contribution shrinks to §4 plus the bounds, and the paper stands or falls on computational scale alone. Plan for that contingency rather than discovering it in review.

**Identity discipline.** The paper should be the E-ansatz paper, not the state-of-the-art-at-size-16 paper; trying to be both bloats it and invites unfavorable comparison with CMSV on their home turf.

**One practical note.** The contribution mix has shifted during review — the ladder, the recursion, the bounds, and the verified tower came from this analysis, not the proposal. Before drafting begins, the authors should settle scope and credit explicitly; it changes what "the basic paper" even is.

Bottom line: 2+3 fused is publishable at the Experimental Mathematics tier, but it is roughly two to four months of disciplined work away — the surjectivity theorem, the literature pass, and the scaled exhaustive run are the three things standing between this outline and a submission, and the first two carry kill-risk for specific claims, not just delay.

# On failure modes

Good — the collaboration context helps; I'll frame these as team workstreams rather than reviewer warnings. And your pushback is partly right on all three, so let me be precise about where each risk actually lives, because in each case it's not quite where the obvious objection points.

## 1. Surjectivity: you're right about the math, but the math isn't the risk

Three statements got bundled under "surjectivity," and they have very different statuses.

**S1: every equal-norm Gaussian quadruple has a multiplicative description.** Trivially true by unique factorization in ℤ[i]. This is your point, and I concede it entirely — no kill risk here.

**S2: the syndrome formalism and W cover all solutions with well-defined invariants.** True, but the write-up must extend the bit-vector syndromes to exponent vectors (repeated primes give eᵢ+1 choices per prime, including the real factor qq̄ = p), handle unit twists and the ramified prime 2, and prove a canonization lemma. Without the lemma the "10 syndrome patterns" census isn't a well-defined statistic — unit and relabeling orbits could merge or split patterns. Provable, a week of careful writing, but it's the fiddly kind of provable.

**S3: the implemented enumerator visits every solution with L₁ ≤ B.** This is the load-bearing statement for every exclusion claim, it's a joint property of math and code, and as described it is currently false. Here's the concrete hole: a norm with a square factor p², p ≡ 1 (mod 4), admits *primitive* solutions — one root pair can use the q² representation while another uses qq̄, so the roots share no common factor and no scaling argument reduces to a squarefree norm. Roughly 39% of integers are non-squarefree, so an "exhaustive to B" claim built on the multiplicity-free algorithm has a hole at a positive density of norms. It's not that the missing solutions are likely to exist — it's that the exclusion sentence "no E₅ with L₁ ≤ B" becomes unprovable as stated, and the failure is silent: the search returns "nothing found" either way.

There's a falsifiable prediction inside this: the parametric families L₁(k) will land on non-squarefree values at infinitely many k, so if the implementation matches the multiplicity-free description, a self-test asking the enumerator to re-find all of [3]'s own family members will fail at those k. I'd run exactly that test first — it either exposes the hole or reveals the code already handles more than the proposal describes.

So: no false mathematics lurking, agreed. The kill target is the headline exclusion claim, and the fix is code plus lemmas, not lemmas alone.

## 2. The literature pass: four questions and a live exhibit

Deliverable, concretely: (a) do the objects exist under another name; (b) are the specific theorems in print; (c) what is [3], with its families in symbolic form; (d) do existing tables overlap the catalog. Each has a specific place it would hide:

**The recursion (Thm 3.1) vs. BMR's litter conditions.** We've only read [1]'s abstract. If their litter conditions already articulate the level-wise structure, Theorem 3.1 shrinks from "new structure theory" to "a congruum-lift reading of [1]" — still useful (the bottom-lift algorithm and MITM duality would remain ours), but the framing changes. Reading the full paper is mandatory before any "we introduce" sentence.

**Classical symmetric parameterizations.** Chernick's 1937 AMM paper is literally titled "Ideal Solutions of the Tarry-Escott Problem," and Gloden's 1944 book is the compendium of multigrade constructions. If [3]'s two families are Chernick's or Gloden's rediscovered — plausible, since equal-pair-sum structure is the natural symmetric ansatz — then parts of §2–3 are 1937/1944 material and the related-work section has to say so before a referee does. Dickson's *History* Vol. 2 and Dorwart–Brown 1937 round out the classical sweep; Caley's thesis covers the ℤ[i] side where our §4 parameterization ideas may be partially anticipated.

**Existing tables.** Chen Shuwen's multigrade site and Meyrignac's eslp database have hosted degree-7 symmetric tables for decades. If a comparable E₄ table already exists there, the catalog's novelty reduces to the syndrome census, exhaustiveness certification, and scale — survivable, but it must be known before writing, not after.

**Foreclosure of direction 5.** CMSV don't just prove required divisors — they construct families of ideal solutions in local settings, which present obstructions to establishing certain arithmetic information about solutions. Translation: some divisibility improvements are *provably impossible* by their local families. Whether the towered subclass escapes those obstructions is exactly what determines if direction 5 is a theorem or a dead end, and only a careful read of their §3 answers it.

**The live exhibit.** During routine searching this session I hit a paper from last month: Tsai et al., motivated by anomaly cancellation in chiral gauge theory, studying the symmetric locus in the ideal degree-three PTE problem, where the symmetry reduces the problem to a sum-of-two-squares equation governed by representations as sums of two squares, complete with an asymptotic count N_sym(H) = (4log2/3π²)H³log H + O(H³). That is the E₃-level of our §8, published, with the same ℤ[i] mechanism, one month old. Two consequences: our density section must start at E₄ and cite them, and — more importantly — it proves this micro-area is active *right now*. The risk isn't only being anticipated by 1937; it's being anticipated by next quarter.

Cost of the pass: two or three focused days on MathSciNet/zbMATH plus reading BMR and Caley in full, plus identifying [3] — which is also a blocking dependency, since R6 (the family-intersection computation) cannot even be specified without the families' symbolic forms.

## 3. Exhaustive search: the FLOPs are engineering; the epistemology isn't

On the compute itself you're right, and the arithmetic is comfortable: within B = 10¹², an admissible L₁ has at most 9 split prime factors (the product of the first ten primes ≡ 1 mod 4 already exceeds 10¹³), so |U| ≤ 2⁸ per norm and the pair count per norm is tiny; the cost is ~10¹⁰–10¹¹ admissible norms, i.e., sieving plus ~10¹⁴ cheap operations. Days on decent hardware. Exhaustive-to-10¹² would beat the unconditional L₁ > 5.8×10⁶ bound by five orders of magnitude — a genuinely respectable headline.

The risk is what "exhaustive" means for a null result. A bug in a *discovery* search produces a false solution, which the verification harness catches in milliseconds. A bug in an *exclusion* search produces nothing, forever — it is invisible in the output by construction. Three specific instances: the S3 hole above (a whole class of norms skipped, silently); every "WLOG" in the pair-ordering and sign conventions, each of which is a small lemma that, if wrong, silently halves the search space; and any ladder-based congruence pre-sieve, whose *completeness* ("we discard these L₁ only because provably no solution lives there") is itself a theorem — an aggressive sieve is exactly how you buy speed by quietly sacrificing exhaustiveness.

The mitigation protocol is standard for record computations and should be in the paper: two independent implementations (ideally two people); synthetic solution injection; the family self-test from §1 including non-squarefree members; cross-checking per-norm counts against the 𝔽ₚ density model, where deviations flag bugs rather than discoveries; pre-registered bounds and checkpointed logs. This is why I resist "just engineering": the deliverable isn't a number, it's a *certificate*, and the certificate's validity is mathematics.

## Bottom line

None of the three risks is likely to kill the project. Each can kill the headline sentence of one section: S3 kills the exclusion claim, the literature pass can kill the "first/new" framing of §§2–4 and possibly direction 5, and search epistemology kills the credibility of the null result. All three are cheap to de-risk now — roughly: one to two weeks for the S2/S3 lemmas plus the exponent-vector extension of the enumerator, two to three days of focused literature work plus two full paper reads, and a hardening protocol that's mostly process. Since you're working with Roi and Gil directly, these partition cleanly into three parallel workstreams, and I'd sequence the literature pass first — it's the cheapest and the only one whose outcome changes what the other two should claim.
