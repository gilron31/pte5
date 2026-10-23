# The algebraic-geometry branch: report and audit

2026-10-23. A self-contained account of the research branch begun on 2026-10-21:
what was done, what came out, what an audit of the work changed, what did not
work, and which leads look worth pulling next. The evidence is in
`ag_research_notes.md` (entries N1–N19), the scripts in `ag_research/`, and the
background in `algebraic_geometry_plan.md` (lines G1–G6, with exercises and
reading). Nothing here has been put into `paper.md`.

Status labels, as in the plan:

- **PROVEN** means an exact computation or a theorem with proof.
- **EVIDENCE** means strong support, but not a proof.
- **NEG** means a negative result within the stated scope.
- **OPEN** means not settled.

---

## 1. Summary

The goal was to attack the $E_5$ problem with algebraic geometry and computer
algebra. An $E_5$ is a 16-root solution of the paper's tower equations,
equivalently a rational point of a 4-fold $V_5 \subset \mathbb P^{15}$ off its
degenerate locus. Six results stand after the audit:

1. **Bottom-lift certificates (G5, PROVEN, after the audit's fix).** For each of
   42 seed forms, no $E_5$ has that seed, up to scaling, as its derived $E_4$, for
   any $\mu$. The 42 forms are:
   - every $E_4$ of the $L_1 \le 10^{12}$ catalog, in both 45° twin forms;
   - Bremner's family 1 for $k \in \{-5, \dots, 4\}$ (nondegenerate members) and
     family 2 for $5 \le \lvert k\rvert \le 10$, with roots up to
     $2 \cdot 10^{27}$.

   The method is genus-2 descent, Chabauty–Coleman and a Mordell–Weil sieve on
   one rank-1 quotient of the seed's lift curve, which has genus 98,305.
   Combined with the catalog, it shows that an $E_5$ needs
   $L_1 > 10^{12}$ and $L_2 > 5 \cdot 10^{11} c^2$.
2. **No $E_5$ from Bremner's families (G4, PROVEN).** No two members of his two
   one-parameter families form an $E_5$: not within family 1, not within family
   2, and not one from each. The proof is exact: implicitized image curves in
   invariant space, then their rational singular points.
3. **No local obstruction (G2, PROVEN).** Congruences alone can never rule out
   $E_5$. CMSV's divisibility floors are locally optimal at $p = 5, 13, 17, 37$,
   and nothing at all is forced at 47.
4. **The geometry of Bremner's slice (G3, PROVEN).** His hyperplane section
   $A - C = E - G$ is the unique plane of its kind. It is a
   $(\mathbb Z/2)^4$-cover of $\mathbb P^2$ branched over four conics, with four
   K3 quotients. Family 1 is the slice cut by one further hyperplane.
5. **A sporadic $E_4$ (G3, EVIDENCE).** The catalog entry
   $L_1 = 33065996717$ lies on the slice but on neither family, and the slice has
   no other new points up to height $10^4$.
6. **$V_4$ is of general type (PROVEN); $V_5$ almost certainly is (EVIDENCE).**
   $V_5$'s worst singularities reproduce the cone over $V_4$, which is the tower
   repeating itself in the local geometry.

**What the audit changed** (§5):

- It found a soundness gap in the Mordell–Weil sieve, which left 9 of the 49
  certificates incomplete as proofs, though none was wrong. After the fix, every
  seed form is re-certified by two independent quotients.
- It corrected an $L_2$ bound by a factor of 2.
- It found that $V_5$'s singular locus is larger than recorded: at least 77,824
  non-canonical points instead of 128, plus a missed locus at $\sigma = 0$. The
  general-type margin falls from about $2 \cdot 10^4$ to about 34.
- It narrowed the scope of the family theorem: family 2's 45° twins are not
  covered.

## 2. The setting, briefly

An $E_n$ is a set of $2^{n-1}$ positive integers, the roots, arranged in a
binary tree whose levels satisfy the tower conditions. Equivalently,
$\prod (x^2 - r_k^2) = E_n(x)$, where $E_1 = x^2 - L_1$ and
$E_{k+1} = E_k^2 - L_{k+1}$. Projectively these are the rational points of a
variety $V_n$ of dimension $n - 1$:

| variety | ambient space | type |
| --- | --- | --- |
| $V_4$ | $\mathbb P^7$ | complete intersection of type $(2,2,2,4)$ |
| $V_5$ | $\mathbb P^{15}$ | complete intersection of type $(2^7, 4^3, 8)$, $K = \mathcal O(18)$ |

There are two ways to build an $E_5$:

- **Top-down.** Two $E_4$ sharing $L_1, L_2, L_3$ glue into an $E_5$.
- **Bottom-up.** An $E_5$ lies above its *derived* $E_4$, with roots
  $\rho_k \propto (r_{2k}^2 - r_{2k+1}^2)/2$ (paper Prop. 3.1, Cor. 3.2). Above
  a fixed seed $\rho$, an $E_5$ is a rational $t = L/\mu$ that puts all 16
  numbers $t \pm \rho_k$ in one square class. That is a rational point on a curve
  $D$ of genus 98,305. Every even-size subset $S$ of the signed roots gives a
  quotient $C_S: y^2 = \prod_{e \in S}(t + e)$, which is independent of $\mu$. So
  one quotient whose rational points are all known settles the seed for every
  $\mu$.

## 3. What was done

| line | question | what was done | result | status |
| --- | --- | --- | --- | --- |
| G1 | what kind of variety is $V_4$? | singular locus | 864 ODPs, all degenerate; $K = \mathcal O(2)$, so general type | PROVEN |
| G1 | and $V_5$? | singular strata (N12), local types (N12a), transversal slices (B1) | $cA_1$ along curves; cone-over-$V_4$ points with $a = -3$; general type by a plurigenus count | EVIDENCE (see §5.3) |
| G2 | can congruences decide $E_5$? | Hensel towers; explicit $\mathbb Z_p$-witnesses | no obstruction at any depth; §4 table settled at 5, 11, 13; CMSV optimal at 5, 13, 17, 37, 47 | PROVEN |
| G3 | where do $V_4$'s rational points live? | catalog to $10^{12}$, identities, special planes, slice model, curve searches, slice point search to height $10^4$ (CUDA) | Bremner's plane is unique; families explained; the exemplar is off-family and isolated | PROVEN / EVIDENCE |
| G4 | can two family members form an $E_5$? | images in $(J_1, J_2)$, exact for family 1, multi-modular plus exact check for family 2; $\varphi$ on the slice | no, in all three combinations; $\deg\varphi = 32 = 8 + 24$, no slice–slice partners up to 400 | PROVEN / EVIDENCE |
| G5 | which $E_4$ carry an $E_5$ above them? | genus-2 descent, Chabauty–Coleman (Sage), MW sieve, pipeline over 42 seed forms | none of the 42 | PROVEN |

Infrastructure along the way:

- a C++/OpenMP and a CUDA slice sieve, the latter 90 s at height 5000 on an
  RTX 4060;
- multi-modular reconstruction for eliminations that stall in characteristic 0;
- a parallel certificate pipeline.

## 4. Notable finds

### 4.1 Per-seed certificates, and the rank signature (G5)

Chabauty–Coleman needs a quotient whose Jacobian rank is below its genus. The
elliptic quotients never qualify: two points of infinite order are built into
the problem, the trivial lift $\mu = 0$ and the seed's own $L_1$ (plan §7.2).
Genus-2 quotients do qualify.

For Symes' $E_4$, 77 of 400 random non-symmetric 6-subsets (19%) have descent
rank bound exactly 1 (N13). Every one of the 42 seed forms was certified by the
*first* rank-1 quotient tried. The cost is a few seconds to five minutes
per seed, and it does not depend on the size of the roots, since $10^{27}$ was
handled.

The converse is the striking part. The derived $E_3$ of Symes' $E_4$ genuinely
lifts, and **every** 6-subset quotient of its lift curve has rank bound at
least 3 (N16): the genuine lift supplies extra independent points. So the rank
distribution of the genus-2 quotients is a cheap signature of liftability, and
it can serve as a screen (lead L3).

The certified statement is about the derived $E_4$, which is a different region
from the classical bound on $L_1$. With the catalog's collision check
(Prop. 6.1), an $E_5$ needs $L_1(E_5) > 10^{12}$ **and**
$L_2(E_5) > 5 \cdot 10^{11} c^2$, where $c$ is the content of the derived
sequence. The audit checked the corrected constant (§5.2).

### 4.2 Bremner's families cannot combine (G4)

The two families are curves on the slice: in its $(s:t:u)$-model, family 1 is a
sextic and family 2 a septic, both of genus 1. An $E_5$ from two members needs
equal projective invariants $(J_1, J_2) = (8L_2/(2L_1)^2,\ 8L_3/(2L_1)^4)$ with
different $L_4$.

- **Family 1:** the image curve was computed exactly, with degree 12 in $J_1$.
- **Family 2:** char-0 elimination stalled for 16 hours. Instead, 11 primes near
  $2^{31}$, CRT and rational reconstruction gave the image: 223 terms with
  145-bit coefficients. It was then proved exact by evaluation in number fields
  at 449 points, above a degree bound of 443.
- Since the maps to the images are 4:1 and 2:1 (each fiber one $E_4$), collisions
  can only sit over singular points of the images. Their rational singular points
  were enumerated, and every preimage is either the same $E_4$, degenerate, or
  irrational.
- The cross-family case is a resultant of degree 332, with two independent
  routes agreeing.

### 4.3 The exemplar (G3)

The catalog $E_4$ with $L_1 = 33065996717$ satisfies Bremner's relation (twin
form) and lies on neither family curve; that is exact, checked over all 32
preimages and all symmetries. Two searches for something that would explain it
came back empty:

- **Linear relations.** LLL finds no second short linear relation, so the
  exemplar is not on a family-1-type hyperplane section.
- **Companion points.** A height search of the slice to $10^4$ finds only the
  two families and the exemplar, at heights 355 and 515.

The working hypothesis is a **sporadic rational point on a surface of general
type**, which Bombieri–Lang allows. Only a new curve of genus at most 1 through
it could overturn that (leads L6, L7).

### 4.4 The structure of Bremner's slice (G3)

In the $(s:t:u)$-coordinates found here, the four square conditions of the
slice become
$$B^2 = 5s^2 + 2st + u^2,\quad D^2 = 5s^2 - 2st + u^2,\quad F^2 = 5s^2 + 2su + t^2,\quad H^2 = 5s^2 - 2su + t^2.$$
Consequences:

- **The cover.** The slice is a $(\mathbb Z/2)^4$-cover of $\mathbb P^2$ over
  four conics. Its 15 intermediate double covers are 4 quadrics, 6 del Pezzo
  surfaces of degree 2, 4 K3 surfaces and one octic double plane of general
  type.
- **Uniqueness.** Among planes of the form $G = aA + bC + cE$, Bremner's is the
  only nondegenerate one on which all four conditions become conics (N6).
- **Family 1** is literally the slice cut by Bremner's second relation
  $A + B + C + D = E + F + G + H$ (N11).
- **No cheap curves.** There is no further genus-0 or genus-1 curve over a line or
  a conic (N8b), so any third family lies over a curve of degree at least 3,
  like the known two (degrees 6 and 7).

### 4.5 $V_5$'s singularities repeat the tower (G1)

At the points where all eight branch divisors of the Kummer cover meet, the
tangent cone of $V_5$ is *exactly* the affine cone over $V_4$, with the same
ideal and degree 32. The audit found two more families of points of this type:
over the vertex of the base cone, and at $\sigma = 0$.

Each such point has discrepancy $a = -3$, so it is not log canonical. Elsewhere
the singularities are transversal $A_1$, which is canonical. The local geometry
reproduces the global recursion one level down.

### 4.6 Local theory (G2)

A Hensel construction gives $\mathbb Z_p$-towers of every depth for every $p$,
so no congruence argument can ever exclude $E_5$. Explicit witness towers built
from catalog seeds show that CMSV's divisibility floors cannot be improved
locally at $p = 5, 13, 17, 37$, and that nothing is forced at 47. The primes
where a better local floor could still exist are $\{2, 3, 7, 11, 19, 23, 29, 31,
41, 43, 53\}$ (N7).

## 5. Audit

I re-read every notebook entry against the scripts that produced it, and
re-derived the arguments behind each PROVEN label. Where a claim rested on code,
I re-ran it, or ran an independent check. Five findings changed something; the
rest held.

### 5.1 A soundness gap in the Mordell–Weil sieve (G5): fixed, all re-certified

**The gap.** The sieve writes a rational point's class as $cD_0 + T$ with
$D_0 = [P_0 - \infty]$. This is valid modulo a sieve prime $q$ only if the index
$m$ of $\langle D_0\rangle + J[2]$ in $J(\mathbb Q)$ is prime to $\#J(\mathbb F_q)$.

- **What was checked:** the code proved $\ell \nmid m$ only for the odd primes
  $\ell$ dividing $\mathrm{ord}(D_0 \bmod q)$.
- **What was needed:** every odd $\ell$ dividing $\#J(\mathbb F_q)$. If
  $\ell \mid m$ and $\ell \mid \#J(\mathbb F_q)$ but $\ell \nmid \mathrm{ord}(D_0)$,
  the sieve can wrongly discard a genuine point.
- **Where it was:** both versions, `mw_sieve.py` (N15) and `lift_pipeline.py`
  (N17). N15's docstring even states the correct condition.

**The fix.** At each sieve prime, multiply all classes by
$k_q = \prod \ell^{v_\ell(\#J(\mathbb F_q))}$, taken over the odd $\ell$ not
dividing $\mathrm{ord}(D_0)$. This projects away the unknown $\ell$-parts.
Saturation is still proved at the $\ell$ dividing the orders.

Why it is sound: let $\Lambda = \mathrm{lcm}_q(\#J(\mathbb F_q)/k_q)$.

- $\gcd(m, \Lambda) = 1$, since $m$ is odd and every odd prime of $\Lambda$ is
  proven saturated.
- So $k_q[P - \infty] = k_q(cD_0 + T)$ holds with one integer
  $c = n m^{-1} \bmod \Lambda$ and one $T$, the same at every $q$.

The projected sieve is weaker, so N16's omission test still passes: survivors
only grow.

**The re-run** (`audit_certs.py`, `audit_certs.log`). All 42 logged certificates
and the 7 N15 quotients were redone:

| check | result |
| --- | --- |
| seeds are nondegenerate $E_4$ (polynomial identity; a perturbed seed fails) | 42 / 42 |
| rank 1 re-proved; Chabauty re-run, same $p_0$ as logged | 49 / 49 |
| old certificates whose sieve primes had untested $\ell$-parts | 16 / 49 |
| of those, relying on an $\ell$ never proven saturated (the old proof was incomplete) | 9 |
| same quotient re-certified by the fixed sieve (default pool; the other 3 with pool to 220) | 46 / 49, then **49 / 49** |
| **independent second certificate** (fresh random quotient, fixed pipeline), per seed form | **42 / 42** |
| point search $H = 3000$ on every quotient: points beyond the known ones | none |

The gap was real. In 9 certificates the old proof rested on an $\ell$-saturation
that was never established; the unproven $\ell$ were 3, 5, 7, 19 and 47. No
conclusion turned out false: after the fix, nothing new survives anywhere.

Three quotients were not re-certified by the fixed sieve with the default prime
pool:

- cat12781559905 (point form), which had relied on an unproven 5-saturation;
- F1k4 (twin form);
- N15's curve [1,9,11,12,14,15], which needed a pool to 220 even originally.

All three re-certify with the pool widened to primes below 220
(`ag_research/audit_retry_wide.py`, log `retry3_wide_pool.log`, 30–112 s each).

**All 49 quotients are certified by the fixed sieve, and every one of the 42 seed
forms has two independent certificates.** N15's theorem for Symes' $E_4$ stands
with all 7 quotients.

### 5.2 The $L_2$ bound was overstated by a factor of 2

The catalog's $L_1$ is the paper's $L_1 = (r_{2k}^2 + r_{2k+1}^2)/2$ of the
*root* form. The point form has odd pair sums and $L_1 = L_1^{\rm cat}/2$.

- **With primitive roots** and $c$ the content of the derived sequence, a
  derived $E_4$ in point form is covered only to $L_1' \le 5 \cdot 10^{11}$. The
  corollary must then read $L_2(E_5) > 5 \cdot 10^{11} c^2$, not $10^{12} c^2$.
  Corrected in N17.
- **In the integral normalization** (point forms doubled so that $L_1 \in \mathbb Z$),
  the original "$L_1' > 10^{12}$" is right.

The difference is purely one of convention, but the published sentence has to
fix one.

Prop. 6.1 is unaffected in the integral normalization: an $E_5$ with
$L_1 \le 10^{12}$ has both halves in the catalog. A new check also closes a
loophole that the campaign's distinct-$L_1$ test left open: the 30 forms have
pairwise distinct $(J_1, J_2)$, so no root form collides with a point form
either.

The paper's own 2-power bookkeeping is still a TODO (§3), and the final
statements should be phrased in its normalization.

### 5.3 $V_5$'s non-canonical locus is bigger than recorded (G1)

N12/N12a counted 128 cone-over-$V_4$ points and skipped $\sigma = 0$ ($L_1 = 0$,
non-real) as degenerate. Degenerate is irrelevant for the Kodaira dimension,
which sees all complex points. Checked by Jacobian ranks:

| locus | count | Jacobian rank | type |
| --- | --- | --- | --- |
| over $d = \sigma(\pm 1, \dots, \pm 1)$ (all 256 sign patterns, not 1) | $256 \times 128 = 32{,}768$ | 7 | cone over $V_4$, $a = -3$ (N12a) |
| over the vertex of the base cone: all roots $\pm 1$ | $2^{15} = 32{,}768$ | 7 | the same (étale over the cone's vertex) |
| $\sigma = 0$, two tree groups of pairs entirely zero | $\ge 6 \times 16 \times 128 = 12{,}288$ | 7 | the same: the local equations are $V_4$'s (one pattern checked at $p = 97$, the other five by the same local computation) |
| $\sigma = 0$, two pairs entirely zero | curves (e.g. over a Fermat quartic) | 10 | rank-4 quadric: $cA_1$, canonical |
| $\sigma = 0$, every pair nonzero | — | 11 (5 of 6 random points) | smooth, except over nodes of $V_4$: the sixth point lies over $d = (c, \dots, c)$, where $V_4$'s Jacobian rank is 3 |

The plurigenus count becomes
$\ge 77{,}824 \times 108\, m^4 \approx 8.4 \cdot 10^6\, m^4$ of conditions,
against $h^0(\mathcal O(18m)) \approx 2.9 \cdot 10^8\, m^4$. That is a margin of
about 34, down from about $2 \cdot 10^4$. **The general-type conclusion survives**,
but the remaining special points (where strata meet, the $K = (0)$ points, the
$\sigma = 0$ points with other zero patterns) now need real discrepancy numbers
before it is a proof (lead L8).

B1 itself held up. Re-run at 6 random points per curve stratum, all 30 give a
transversal ODP. Since the length of the singular scheme is upper
semicontinuous, the $\mathbb F_p$ result implies the char-0 one, which makes B1
a proof once that argument is written down.

### 5.4 The family theorem covers family 2 on the slice, not its twins (G4)

The G4 theorem is about the rational points of the two family curves in the
$(s:t:u)$-plane. Family 1 is closed under the 45° twin, so it is complete.
Family 2's twins are **not on the slice at all** (checked for $k = 5$–10), so
an $E_5$ involving one is not excluded by G4. It is excluded by G5 only for
$k \le 10$, and only in the derived-$E_4$ role.

The theorem's inputs did check out: every family member tested lands on its own
curve (4 or 2 of its 8 signed labelings, matching the 4:1 and 2:1 maps) and
never on the other.

### 5.5 What the results rest on

These are the assumptions that remain. None was found wrong.

- **Catalog completeness to $10^{12}$:** the S-A enumeration in `code/native`.
  This branch only checked that each entry is an $E_4$.
- **Paper Prop. 3.1 and Cor. 3.2:** the bottom-lift criterion. Its use here was
  validated: Symes' known lift lies on every quotient of its $E_3$-level lift
  curve (N16).
- **Software:**
  - Sage's Coleman integrals, which passed the Weierstrass-point identity to
    $O(17^{15})$;
  - PARI's `hyperellcharpoly` and `factor`;
  - Singular's Gröbner bases, over $\mathbb Q$ and $\mathbb F_p$.
- **Descent targets:** the local image sizes $2^4$, $2^6$, $2^2$ (standard for
  full rational 2-torsion). The certifier refuses any curve whose sampled image
  falls short, so a short sample costs completeness, not soundness.
- **Chabauty tail control:** the Strassmann bound reads coefficients to $t^{40}$,
  and the tail beyond is bounded by the integrality of the local expansions at a
  good prime. That integrality is standard but was not re-verified separately
  for Sage's local coordinates.
- **Mod-$p$ computations standing in for char 0:** $\varphi$'s degree (5 fibers),
  B1, and the slice sieve. The exact parts, the G4 images and the G5
  certificates, do not depend on this.

## 6. Directions that did not work

| approach | outcome | why, and the lesson |
| --- | --- | --- |
| elliptic quotients of rank 0 (the cheapest G5 certificate) | never: every elliptic quotient has rank $\ge 1$ | two infinite-order points are built in ($\mu = 0$ and $X = 2L_1$). Structural, so no seed reopens it; go to genus 2 |
| calibrating G5 on a real lift ($E_3 \to E_4$) | impossible: all quotients have rank $\ge 3$ | genuine lifts add independent points. Replaced by the omission test; it became the rank signature (§4.1) |
| deciding $E_5$ by local / $p$-adic arguments | provably impossible | $\mathbb Z_p$-towers exist at every depth. Only lower bounds at 11 specific primes can still improve |
| integer seeds as $p$-adic witnesses at mid-size $p$ | weak ceilings | needs all 16 values to be squares mod $p$ ($2^{-16}$); use $\mathbb Z_p$-seeds (E1) |
| lines and conics on the slice (genus $\le 1$ curves) | none rational and nondegenerate | any third family needs degree $\ge 3$ in $(s:t:u)$ |
| family-2-type cubic through the exemplar | none (6 split lines, coefficients $\le 15$, 32 preimages) | narrow ansatz; the exemplar only sits on genus-1 curves of each K3 separately |
| LLL for a second linear relation at the exemplar | none | rules out family-1-type (hyperplane) curves through it |
| height search for companions of the exemplar | nothing new to $10^4$ | supports "sporadic", which cannot be proved this way |
| $\varphi$-scans for slice–slice $E_5$ | no partner in 240,000 base points | consistent with Hilbert irreducibility; needs the geometry of $X'$ (D1) |
| char-0 elimination of family 2's image | stalled for 16 h | multi-modular with an exact check took 51 s |
| primary decomposition for $V_5$'s local types | timed out after an hour | random $\mathbb F_p$-points plus transversal slices took 5 s |
| Singular `modStd`; a naive MW-sieve loop; a sieve pool to 400 | memory blow-up; hang; too slow | sequential modular; incremental greedy sieve; pool $\le 120$ |
| the first MW-sieve implementation | unsound in a corner case (§5.1) | write the soundness condition into the code, then test the code against it |

## 7. Leads for future work

Ranked by expected value over cost.

- **L1. Finish G5's rigor (cheap).**
  - Refine the projection to kill only the $\ell$ that cannot be proven saturated,
    which would re-certify the one quotient lost in §5.1 with more sieve power.
  - Add degree-2 divisors to the 2-adic sampling (A4).
  - Have one or two certificates re-checked by an independent implementation
    (Magma's `Chabauty` and `MordellWeilSieve`, if access appears).
- **L2. Scale G5 (cheap per seed).** Both families grow quadratically in digits.
  Family 1 reaches 57 digits at $k = 6$ and family 2 reaches 40 digits at
  $k = 20$, still within PARI's factoring of the root differences, which the
  descent needs. A catalog to $10^{13}$–$10^{14}$ from the S-A campaign could be
  certified at minutes per seed; the bottleneck is generating the catalog.
  - Beyond about 80 digits, the differences must be factored algebraically
    through the family parametrization.
- **L3. The rank signature as a liftability screen.** The descent takes 0.04 s
  per curve. A seed with a genuine lift has no rank-1 genus-2 quotient (§4.1). So
  surveying about 20 quotients per seed sorts a large catalog into "certifiable
  now" and "suspicious". Any seed that resists certification is the interesting
  one.
- **L4. Top-side certificates.** G5 excludes an $E_4$ as the *derived* $E_4$;
  the other role, *half* of an $E_5$, is open beyond the catalog. Partners of a
  fixed $E_4$ are rational points on the genus-21 curve
  $E_3(x) + E_3(x') = 0$ (plan §3.3), which has many involutions. Look for a
  quotient with rank below its genus. A partner would add 128 new points. The
  seed's own roots already supply 128 known points ($E_3(x) = \pm s$), which may
  force the ranks up, so compute rank bounds first.
- **L5. Close G4's scope.** Compute the image of twin(family 2) in
  $(J_1, J_2)$ with `modular_image.py` and intersect it with both family images.
  It is the same machinery as N8d.
- **L6. K3 fibrations through the exemplar (C2).**
  - Take the elliptic fibrations of the four K3 quotients through its points.
  - Compute the fiber ranks.
  - Intersect the candidate curves.

  This is the cheapest test of "sporadic".
- **L7. Picard numbers of the slice's K3s (C3).** Point counts over
  $\mathbb F_{p^n}$, now feasible on the GPU. They give the lattice of curve
  classes and a finite list of low-degree candidates for a third family.
- **L8. Make "$V_5$ is of general type" a theorem (B2).**
  - Enumerate the full singular locus, including $\sigma = 0$ and all zero
    patterns.
  - Compute discrepancies at the special points.
  - Redo the count: the margin is about 34, so it is not automatic.
- **L9. The degree-24 correspondence $X'$ (D1).** Its genus over lines and its
  Kodaira dimension decide whether slice–slice $E_5$ candidates can come in
  families.
- **L10. Local floors (E1).** $\mathbb Z_p$-seed witnesses at
  $2, 3, 7, 11$ and the mid-size primes, the only places a better local bound can
  still come from.

## 8. Where things are

- **Notebook:** `ag_research_notes.md`; N19 is this audit.
- **Plan and status table:** `ag_research_plan.md`.
- **Guide:** `algebraic_geometry_plan.md`.
- **Scripts:** `ag_research/`; the index is in `ag_research/README.md`.
- **Audit run:** `ag_research/audit_certs.py` (log `audit_certs.log`).
- **Fixed pipeline:** `ag_research/lift_pipeline.py`, where the `mw_sieve`
  projection is on by default.
