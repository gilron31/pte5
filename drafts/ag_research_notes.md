# AG research notebook

Working notes for the research branch of `algebraic_geometry_plan.md` (lines
G1–G6). Not paper text. Each entry: what was run, what came out, what it means,
what next. Scripts are in `drafts/ag_research/`; every number below was produced
by them.

Notation: Bremner's labels, pairs $(A,B),(C,D)\mid(E,F),(G,H)$. The **slice**
$S_B = V_4 \cap \{A - C = E - G\}$. A catalog entry has two 45°-twin forms: the
*root form* (roots $m \pm n$) and the *point form* (roots $m, n$).

---

## 2026-10-21

### N1. The catalog to $10^{12}$ (S-A rerun, scratchpad + `catalog_1e12.csv`)

`enum_campaign --bound=1e12`: 17.6 billion admissible norms, 725 meet
candidates, 710 imprimitive, **15 $E_4$**, all certified, 10/10 known recovered,
558 s on 16 threads. New beyond $10^{11}$: $L_1 = 127393252597$,
$699758104565$, $806692500965$. All 15 have distinct $L_1$, so there is no
top-collision $E_5$ with $L_1 \le 10^{12}$ (Prop. 6.1, now to $10^{12}$).

### N2. Who is who among the known $E_4$ (`known_e4.py`, `census.py`)

Pooled: the catalog, Bremner's Table 1 (11 rows) and both families at
$|k| \le 12$, each in both twin forms. Identities found:

| catalog $L_1$ | identity |
| --- | --- |
| 67405 | Table 1 row 1\* = family 2, $k = 5$ (twin) |
| 600445 | row 2\* = family 1, $k = 1$ (and $k = -2$, twin) |
| 5915065 | row 3\* = family 2, $k = \pm 6$ |
| 43937864165 | row 9\* = family 2, $k = \pm 7$ (twin) |
| 55445869, 81524605, 166045477, 242048509, 532946005 | rows 4, 5, 6, 7, 8 (Bremner's non-starred) |
| 12781559905, 17247059713, **33065996717**, 127393252597, 699758104565, 806692500965 | not in Bremner's tables or families |

Rows 10\* and 11\* are family 2 at $k = \pm 8, \pm 9$, so *every* starred row of
Table 1 is a family member.

Two structural facts about the families:

- **Family 1 is closed under the 45° twin:** $\mathrm{twin}(F_1(k)) = F_1(-k-1)$
  for every $k$ tested.
- **Family 2 has $F_2(k) = F_2(-k)$** (Bremner's remark), and its twins are not
  family members.

### N3. Small linear relations (`census.py`)

Test every signed zero sum $\sum c_i \lvert r_i\rvert = 0$, $c \in \{-1,0,1\}^8$,
at least 3 terms, on both forms. Recorded as terms per pair.

- A (1,1,1,1) relation (one root from each pair) is, up to the signed tree
  symmetries, exactly Bremner's $A - C = E - G$. A (2,2,2,2) relation in one form
  is the same relation in the twin form, since $(A+B) \pm (A-B) \in \{2A, 2B\}$.
  So there is **one** relation class, "on the slice", and both families live on
  it. Family 1 is on it in both forms, which is Bremner's (2-3) and (2-4).
- **On the slice:** the 4 family members, and **the exemplar
  $L_1 = 33065996717$ (twin form)**.
- **Off every ±1 relation:** the 5 non-starred Bremner rows and the other 5 new
  catalog entries.
- Small entries (67405, 600445) carry extra accidental relations; ignore them.

So 5 of the 15 $E_4$ with $L_1 \le 10^{12}$ lie on the slice. For a generic $E_4$
with roots near $10^5$ the chance of a ±1 four-term relation is about $10^{-3}$,
so the exemplar on the slice is very unlikely to be an accident.

### N4. The exemplar is on the slice but on neither family curve (`invert.py`)

Inverting Bremner's parameterization ($(\lambda:\mu:\nu:\phi) = (s:B:u:D)$, then
the quadric $\mu^2 + \phi^2 = 2\nu^2 + 10\lambda^2$ back to $(x:y:z)$) over all
signed labelings gives **32 preimages** of the exemplar in the $(x:y:z)$-plane.
The smallest are $(\pm 96, \pm 355, 225)$ and $(\pm 110, \pm 515, 297)$.
None lies on FAM1, FAM2, Bremner's rank-0 cubic or his degenerate lines.

Control: family member 43937864165 has 16 preimages, and exactly those that are
Bremner's own points lie on FAM2, e.g. $(-177, -201, 11)$, which matches his
Table 3 at $k = 7$. The preimage set is closed under the slice's symmetries, so
"no preimage on FAM1 or FAM2" means the exemplar is on no image of the families
under the known symmetry group either.

**Open question Q1.** Is the exemplar on a third special curve of $S_B$ (a new
family, or a rank-0 curve), or is it a sporadic point of a surface of general
type?

**How family 2 sits in the double planes** (sympy):
$$q_1 - F_2^2 = -(y+z)\,\big((x-y-2z)^2 - 8z^2\big)\,\mathrm{FAM2},\qquad
q_2 - H_2^2 = 4z(x-z)(y+z)\,\mathrm{FAM2},$$
where $F_2, H_2$ are Bremner's explicit cubic forms. So FAM2 is a **common
component of the two sextics $q_1 = F_2^2$ and $q_2 = H_2^2$**. The other
components are degenerate split lines and a line pair defined over
$\mathbb Q(\sqrt 2)$.

That suggests an **ansatz for Q1.** Look for cubic forms
$\Phi = \pm F_y + (y+z)Q$ and $\Psi = \pm H_y + (y+z)Q'$, where $F_y, H_y$ are
the square roots of $q_1, q_2$ on the split line $y = -z$ and $Q, Q'$ are
quadrics. Impose $\Phi(P_0) = \pm F(P_0)$ and $\Psi(P_0) = \pm H(P_0)$ at the
exemplar's preimage $P_0 = (-96, -355, 225)$, which is one linear equation on
each 6-dimensional $Q$, $Q'$. Then test whether $(q_1 - \Phi^2)/(y+z)$ and
$(q_2 - \Psi^2)/(y+z)$ share a cubic factor. The short candidates come from LLL.
Family 2 itself is one solution, and that is the calibration. Not yet run.
(Family 1 is *not* of this type: its $F$ is the rational function
$(T^2 + q_1 - q_2)/2T$.)

### N5. The slice up to height 500 (`slice_search.py`, `slice_search_big.py`)

Exact search of $q_1 = \square$, $q_2 = \square$ over primitive $(x:y:z)$,
$\lvert x\rvert, \lvert y\rvert, z \le 500$. The int64 evaluation is exact since
$240 \cdot 500^6 < 2^{63}$. It takes 21 s.

- **762,144 points; 6 nondegenerate $E_4$:** five family members (67405, 600445
  in both twin forms, 5915065, 43937864165) and the exemplar at $(-96, -355, 225)$.
  **Nothing new.** The exemplar is isolated at this height.
- **The degenerate points lie on split curves.** Both sextics restrict to squares
  along:
  - the **six lines** $z = 0$, $x + y = 0$, $x \pm z = 0$, $y \pm z = 0$
    (Bremner's four degenerate factors plus their mirrors);
  - the **four conics** $x^2 - xy + xz + yz - 10z^2$,
    $x^2 - xy - xz - yz - 10z^2$, $y^2 - xy + xz + yz - 10z^2$ and
    $y^2 - xy - xz - yz - 10z^2$ (signs normalized; `degen_curves.py`), on which
    the two halves of the tower coincide ($L_4 = 0$).

  These are classes in $\mathrm{NS}(X_1)$ and $\mathrm{NS}(X_2)$, input for G3(a).

### N6. Bremner's plane is the only one (`special_planes.py`)

Solved for every plane $G = aA + bC + cE$ in the ACEG model on which the
denominator $A^2 + C^2 - E^2 - G^2$ divides the numerator of $b = B^2$. The
real solutions are:

- the 8 planes $G = \pm A \pm C \pm E$, which are Bremner's slice up to sign
  symmetry;
- $G = \pm A$ and $G = \pm C$, which force a repeated $\lvert\text{root}\rvert$
  and so are degenerate.

There are also two complex solutions. Planes outside this chart reduce to these
by the half-swap. **So up to symmetry Bremner's slice is the unique
nondegenerate plane** on which all four square conditions become conics. G3(e)
is closed: no other "Bremner-type" slices exist.

**The slice in its own coordinates.** On the plane, with $(s:t:u)$,
$$B^2 = 5s^2 + 2st + u^2,\ D^2 = 5s^2 - 2st + u^2,\ F^2 = 5s^2 + 2su + t^2,\ H^2 = 5s^2 - 2su + t^2.$$
So $S_B$ is the $(\mathbb Z/2)^4$-cover of $\mathbb P^2$ branched along four conics.

- The conic pairs $(b, d)$ and $(f, h)$ are each tangent at one point:
  $b - d = 4st$, and on $s = 0$ both are $u^2$.
- The other pairs meet in four non-real points; for example
  $b - f = (t - u)(2s - t - u)$.
- The 15 intermediate double covers are: 4 quadric surfaces (one conic each);
  6 quartic double planes (two conics: del Pezzo of degree 2, rational);
  **4 K3 surfaces** $w^2 = bdf$, $bdh$, $bfh$, $dfh$ (three conics, a sextic with
  simple singularities); and 1 octic double plane $w^2 = bdfh$, of general type.
- Symmetries: $t \leftrightarrow u$ swaps $(b, d) \leftrightarrow (f, h)$, and
  $s \mapsto -s$ swaps $b \leftrightarrow d$ and $f \leftrightarrow h$.

This is cleaner than the $(x:y:z)$ picture of the plan, where two conics have
already been parameterized away.

### N7. How far local arguments can go at the top node (`top_node_witness.py`)

For each prime $p \le 53$ I searched the 15 catalog seeds (both forms; integral
form only at $p = 2$), with $v_p(\mu) \le 4$ and $L < 200$. The output is the
smallest $v_p(\sqrt{L_5})$ over the explicit $\mathbb Z_p$-towers found.
A tower with $v_p(\sqrt{L_5}) = t$ proves that no $p$-adic argument can force
$p^{t+1}$.

| $p$ | CMSV $e_p$ | witness ceiling | |
| --- | --- | --- | --- |
| 5, 13, 17, 37 | 4, 2, 1, 1 | 4, 2, 1, 1 | **CMSV is locally optimal** |
| 47 | 0 | 0 | **nothing forced** |
| 2, 3, 7, 11 | 10, 6, 3, 2 | 29, 14, 10, 4 | open |
| 19, 23, 29, 31, 41, 43, 53 | 1, 1, 1, 0, 1, 1, 1 | 2, 2, 8, 1, 8, 2, 8 | open (search too small) |

The weak ceilings at mid-size $p$ are an artifact. With $v_p(\mu) = 0$ the trick
needs all 16 values $L \pm \rho_k$ to be squares mod $p$, which happens with
chance $2^{-16}$ per $(L, \text{seed})$. Better witnesses need $p$-adic seeds,
not integer ones. Any prime with a cap-0 $\mathbb F_p$-tower (the paper's
"distinct roots" model is nonzero at 4 primes $\le 113$) has a Hensel-lifted
tower with all node valuations 0, so nothing at all is forced there.

**Consequence.** Any improvement of the $E_5$ lower bound by a better *local*
top-node divisor can only come from $p \in \{2, 3, 7, 11, 19, 23, 29, 31, 41, 43, 53\}$.

### N8. G4, rigorous family intersection: method found, family 2 pending

The mod-$p$ census of the plan (§6.3) says no family pair forms an $E_5$, but
that is evidence. For a proof: implicitize each family curve's image in
$(J_1, J_2) = (8L_2/(2L_1)^2,\ 8L_3/(2L_1)^4)$, intersect the two images, and
find the rational points.

- **In $(x:y:z)$ it is too heavy.** `implicit_fam2.sing` ran 20+ min; the
  `modStd` variant forked 9 workers at about 2.5 GB each and had to be killed
  (do not run `modStd` unattended here; `setcores` did not cap it).
- **In $(s:t:u)$ it is cheap**, because $L_1, L_2, L_3$ have degree 2, 4, 8
  there, against 6, 12, 24. The family curves in that plane (`fam_stu.sing`,
  0.1 s) are:
  - **family 1:** an irreducible **sextic**, invariant under $t \leftrightarrow u$
    and $s \mapsto -s$;
  - **family 2:** a **septic** with no symmetry.

  Both are genus-1 curves, so heavily singular.
- **Family 1's image** in the $(J_1, J_2)$-plane: one polynomial of degree 12 in
  $J_1$ and 5 in $J_2$ (`implicit_fam1_stu.sing`, under 1 s; saved to
  `fam1_image_J.txt`). Family 2 is running.

### N8a. No $E_5$ from two members of family 1 (`sing_image.sing`, `twin_check.py`)

Family 1's image curve in the $(J_1, J_2)$-plane is an irreducible curve of
degree 12. The map from family 1's sextic onto it is generically **4:1**: the
fiber over the image of members $k = 1, 2, 3$ is 4 rational points, one $E_4$,
which is the curve's own symmetry group. So two family-1 points with equal
$(J_1, J_2)$ and different $L_4$ can only lie over a **singular point** of the
image.

The affine singular scheme has length 22: 6 rational points and conjugate
clusters of sizes 2, 2 and 12. Above the six rational ones:

- **Five** have fibers made of conjugate pairs, so no rational points at all.
- **One**, $(J_1, J_2) = (4849/4352,\ 182761425/2575826944)$, has 4 rational
  preimages $(1 : \pm 15/11 : \pm 35/11)$ (and the $t \leftrightarrow u$ swaps).
  They are a single $E_4$ (same $j_4$), and not a slice point ($b$ is not a
  square).

**Proposition (computational).** No two rational points of family 1's
$(s:t:u)$-curve with $L_1 \ne 0$ share $(L_1 : L_2 : L_3)$ with different $L_4$.
In particular no $E_5$ is made of two family-1 members, for any $k$. Degenerate
loci are excluded throughout: $s = 0$ gives $A = C$, and $L_1 = 0$ is removed
by saturation.

### N8c. No $E_5$ from one member of each family (`cross_family.gp`, 7 min)

Pull family 1's image curve back along $\varphi$. The result has degree 48 in
$(t, u)$ after stripping powers of $L_1$. Intersect it with family 2's septic.
The resultant in $u$ has degree 332 in $t$: 30 distinct factors, the largest of
degree 222, and two linear. The only rational points of family 2's curve whose
$\varphi$-image lies on family 1's image are $(1:0:-2)$ and $(1:-2:0)$. Both are
degenerate ($\lvert A\rvert = \lvert C\rvert = \lvert G\rvert$).

**So no $E_5$ is made of one member from each family.** Cross-check by an
independent route (`cross_images.gp`): intersecting the two image curves
directly (resultant in $J_2$ of degree 160) gives two rational points. The first,
$(32/25, 0)$, lies under the degenerate pair above ($L_3 = L_4 = 0$) on both
curves. The second, $(14446/12675,\ 539392/11567205)$, has only conjugate
preimages on both curves.

### N8d. No $E_5$ from two members of family 2, so G4 is closed

**Family 2's image curve.** The char-0 elimination stalls (16+ h of wall time,
killed), but over $\mathbb F_p$ it takes 1 s. So `modular_image.py` runs it
over primes just below $2^{31}$, CRTs the coefficients and rationally
reconstructs them. Validated on family 1: identical to the exact char-0 result.
For family 2 it stabilizes after 11 primes:
- 223 terms, degree 28 in $J_1$ and 14 in $J_2$, coefficients up to 145 bits,
  irreducible over $\mathbb Q$ (Singular);
- **proved exactly** (`exact_check_fam2.gp`, 33 s). At each $t_0$, work in
  $\mathbb Q[u]/(G_2(t_0, u))$ and evaluate. The remainder of
  $L_1^{112} F_2(J)$ modulo the septic has $t$-degree at most 443: each
  $u^6$-substitution lowers the $u$-degree by at least 1 and raises the total
  degree by at most 1, and the septic is monic in $u$ up to $-25$. Exact
  vanishing at 449 values of $t_0$ therefore makes it identically zero.

**Twins.** The map from family 2's septic onto its image is generically **2:1**:
two rational points forming one $E_4$, checked at $k = 5, 6, 7, 8$. Its rational
singular points come from `sing_rational.gp`: the gcd of
$\mathrm{Res}_{J_2}(F, \partial_{J_2}F)$ and
$\mathrm{Res}_{J_2}(\partial_{J_1}F, \partial_{J_2}F)$, then rational roots,
calibrated to reproduce family 1's six. There are **7**:

- $(2, 0)$: no preimage with $L_1 \ne 0$;
- $(32/25, 0)$ and $(17/25, 9/200)$: rational preimages, but each pair is one
  $E_4$ with $L_4 = 0$, which is degenerate;
- the other four: only conjugate preimages.

**Theorem (computational, exact).** No $E_5$ consists of two $E_4$ from
Bremner's two families, whether both from one family or one from each. Stronger:
no two nondegenerate rational points of the families' $(s:t:u)$-curves have
equal $(L_1 : L_2 : L_3)$ unless they are the same $E_4$. (N8a, N8c, N8d.)
This closes §6.7 of the plan.

### N4a. Q1 ansatz run (`ansatz_search.py`, `ansatz_exemplar.sh`)

- **Calibration.** At family-2 point $(-177, -201, 11)$, the search over integral
  quadrics with coefficients $\le 10$ recovers FAM2, with $q_2$ a square on it.
  It also finds a cubic on which only $q_1$ is a square, which is a curve of the
  K3 surface $X_1$ but not of the slice.
- **The exemplar.** Over all 32 preimages with coefficients $\le 12$: no
  family-2-type cubic. One preimage, $(-110, -515, 297)$, carries a cubic on
  which $q_1$ is a square but $q_2$ is not. 20 of the 32 preimages admit no
  integral $Q$ at all (the divisibility by $y_0 + z_0$ fails).
- **So the narrow ansatz is negative.** Still open: rational $Q$, larger
  coefficients, other split lines, and curves of family-1 type.

### N8b. Nothing of degree $\le 2$ in the $(s:t:u)$-plane (`tangent_lines.py`, `bitangent_conics.sing`)

A plane curve $\Gamma$ along which all four of $b, d, f, h$ restrict to rational
squares, except $k$ of them, carries a curve of $S_B$ of genus $1 + 2^{k-1}(k-2)$
if $\Gamma$ is a line or conic with distinct branch points: genus 0 for $k = 1$,
genus 1 for $k = 2$.

- **Lines tangent to $\ge 2$ conics: 17 over $\bar{\mathbb Q}$.**
  - The rational ones are $s = 0$, tangent to all four: degenerate ($A = C$).
  - Also rational are the four lines $2s \pm t \pm u = 0$, tangent to two conics.
    On each of them the other two conics restrict to the *same* quadratic, so
    $D^2 = H^2$ (or similar): degenerate.
  - The rest are defined over $\mathbb Q(\sqrt 2)$ or $\mathbb Q(\sqrt 5)$.
- **Conics bitangent to three of the conics**, with rational square restrictions
  (equivalently $C = \lambda b + \ell^2$ for each of the three). The system is
  zero-dimensional with 12 components. The only one over $\mathbb Q$ is
  $C = 6s^2 + t^2 + u^2 = 2L_1$, which has no real points (it is the tautology
  $A^2 + B^2 = 2L_1$). The rest are over $\mathbb Q(\sqrt 2, \dots)$.

**So any further genus-0 or genus-1 curve on the slice lies over a curve of degree
$\ge 3$ in the $(s:t:u)$-plane**, like the two known families (degrees 6 and 7).

### N9. The forgetful map restricted to the slice (`phi_degree.py`, `phi_scan*.gp`)

On the slice the $E_4$ is a function of $(s:t:u)$, and $L_1, L_2, L_3$ are
polynomials there. So $\varphi: \mathbb P^2_{stu} \to \mathbb P(1,2,4)$ is a
**generically finite map between surfaces**. An $E_5$ with both halves on the
slice is a pair $P \ne Q$ in one fiber with different $L_4$ (and $L_1$ ratio a
square).

- **$\deg\varphi = 32$** (Singular, fiber through $(1:2:5)$ and four other random
  points). A generic fiber is **8 rational points**, which are one $E_4$ (its
  signed relabelings: $t \mapsto \pm t$, $u \mapsto \pm u$, $t \leftrightarrow u$),
  plus **one irreducible degree-24 component over $\mathbb Q$**.
- **Special fibers** have 16 rational points that are still one $E_4$. They lie
  on the lines $u = \pm 3t \pm 2s$ and their images, where an extra involution
  fixes the $E_4$. Examples: $(1:3:7)$, $(1:5:13)$, and 192 points in the scan box.
- **Scan.** For every base point $(1:t:u)$ with $\lvert t\rvert \le 150$,
  $1 \le u \le 150$ (about 34,000 points), I looked for rational fiber points with
  a different $E_4$, via the resultant in $u$ and its rational roots. This was
  validated against Singular's counts at three points. Result: **none**.
- **The known slice points** (`member_partners.gp`): the exemplar's point
  $(74250 : 135735 : 120953)$ and every family member with $\lvert k\rvert \le 7$.
  Each has exactly its 8 trivial fiber points and **no rational point with a
  different $E_4$**. That includes $(s:t:u)$-heights up to $10^{76}$, where the
  scan box says nothing.
- **Reading.** By Hilbert irreducibility, a rational slice point has a rational
  "twin" only over a thin set. Slice–slice $E_5$s need the degree-24
  correspondence $X'$ to have rational points, which needs either rational curves
  on $X'$ or sporadic points. Next: the geometry of $X'$ (its genus over lines,
  its Kodaira dimension).

### N10. The slice to height 1200 (`slice_sieve.py`)

A QR-sieve on 14 primes, exact checks on survivors, the six degenerate split
lines masked out. At $H = 500$ it was validated against the exact search (same
762,144 points and 6 $E_4$). At **$H = 1200$**: 2,296 points survive and
**8 nondegenerate $E_4$**: the previous six plus family-2 members $k = \pm 8$
(Table 1 row 10\*) and $k = \pm 9$ (row 11\*). **Nothing new.**

Through height 1200 the slice's rational points are the two families and the
exemplar, which appears at heights 355 and 515. Under Lang, a sporadic point
on a general-type surface is allowed; under the "third curve" reading, more
points on that curve must exist. The families' next members appear at roughly
these heights, so the missing companions are a mild point for the sporadic
reading, not a decisive one.

### N11. Q1: the exemplar has no second short linear relation (`relation_lll.py`)

LLL on the integer relation lattice of the 8 signed slice roots (rank 7,
covolume about $N$, so typical minima about $N^{1/7}$):

| point | $N^{1/7}$ | successive minima (norms) |
| --- | --- | --- |
| family 1, $k = 1, 2, 3$ | 2.7, 16, 233 | **2, 2**, then 3–4, 34–55, 1289–2548 |
| family 2, $k = 6, 7$ | 3.0, 5.9 | **2**, then 3–6, 6–10 |
| **exemplar** | 5.9 | **2**, then 6, 6, 6, 9, 12, 13 |

Family 1 has two anomalously short relations: the slice and Bremner's
$A+B+C+D = E+F+G+H$. I checked symbolically that Bremner's $T = F + H$ is
exactly $A + B + C + D - E - G$, so **family 1 is literally slice $\cap$
hyperplane**. Family 2 and the exemplar have only the slice relation. So the
exemplar is *not* on a family-1-type curve (a hyperplane section); any curve
through it is of the nonlinear, family-2 kind.

### N12. G1.2: the singular locus of $V_5$ (`v5_singular.py`, 0.6 s)

**Criterion** (chain rule through $r \mapsto s$). Work at a point over
$(\sigma : d)$ with $\sigma \ne 0$ whose zero leaves are $r_{2k+1} = 0$ for
$k \in K$, so $d_k = \sigma$ there (signs are a symmetry). $V_5$ is singular
there iff the $4 \times (9 - \lvert K\rvert)$ matrix
$$\big[\ N(d)_{\cdot,k\notin K}\ \big|\ \textstyle\sum_{k\in K} N(d)_{\cdot,k}\ \big]$$
has rank below 4, where $N(d)$ is the Jacobian of $V_4$'s four equations.
$\sigma = 0$ means $L_1 = 0$, which is degenerate and non-real, so it is skipped.

Over the 21 orbits of $K$ under the tree group (order 128), in the base
$\mathbb P^8$:

| $K$ | singular stratum |
| --- | --- |
| $\emptyset$ | **864 lines**, the cone over $\mathrm{Sing}(V_4)$ (the expected value) |
| $(0)$ | 16 points |
| $(0,2)$, $(0,4)$, $(0,1,2,3)$, $(0,1,4,5)$, $(0,2,4,6)$ | curves of degree 16, 80, 16, 16, 8 |
| all 8 | 1 point, $d = \sigma(1,\dots,1)$, with 128 points of $V_5$ above it |
| the other 13 | empty |

**So $\mathrm{Sing}(V_5)$ has dimension 1, is degenerate, and equals the 864
cone lines plus five families of curves plus isolated points.** Over generic
points of the cone lines the cover is étale, so the singularity is a 3-fold ODP
times a curve, which is canonical. Still to do for general type: the transversal
type along the $K \ne \emptyset$ curves (check that the singular scheme is
reduced along them), and log discrepancies at the isolated points, chiefly the
128 points where all 8 branch divisors meet. $K = \mathcal O(18)$ leaves a lot of
room, so general type is very likely; this narrows the proof to finitely many
local checks.

### N4b. Q1 ansatz on all six split lines (`ansatz_exemplar_all.sh`)

I generalized the family-2 ansatz to every split line $\ell$, with
$\Phi = \pm R_1 + \ell Q$ and $R_1$ the square root of $q_1$ on $\ell$; all six
roots were verified with sympy.

- **Calibration.** From its own point, FAM2 is recovered through **4 of the 6
  lines**, so the search is robust to which line a curve is anchored to.
- **The exemplar.** Over all 32 preimages × 6 lines with coefficients $\le 15$:
  12 cubics through preimages on which $q_1$ is a square. On none of them is
  $q_2$ a square. **No family-2-type curve.** The exemplar does lie on genus-1
  curves of each K3 surface separately.

### N12a. G1.2: what the singularities of $V_5$ look like (`v5_worst_point.sing`, `v5_k0_local.sing`)

**The 128 worst points are cones over $V_4$.** These sit over
$(1 : 1, \dots, 1)$, where all 8 branch divisors meet. Put $d_k = 1 - 2y_k^2$:
the local equations become $V_4$'s own equations evaluated at the squares. The
**tangent cone is exactly the affine cone over $V_4$** (degrees $(2,2,2,4)$,
dimension 4, degree 32, identical ideal). The Kummer tower repeats one level
down.

- **Discrepancy.** Blowing up gives $E \cong V_4$ with $E|_E = \mathcal O(-1)$ and
  $K_E = \mathcal O(2)$, so $a = -3$: **not log canonical**.
- **Cost.** Sections of $mK$ on a resolution must vanish to order $3m$ at these
  points. That is about
  $128 \cdot 32 \cdot (3m)^4/4! \approx 1.4 \cdot 10^4\, m^4$ conditions, against
  $h^0(V_5, \mathcal O(18m)) \approx 65536 \cdot (18m)^4/4! \approx 2.9 \cdot 10^8\, m^4$.
  These points cannot destroy general type.

**The 16 "$K = (0)$ points" are not isolated.** In $d$-coordinates they are
$d = (1, -a^2, a^3, a, 0,0,0,0)$, $a^4 = -1$, and its conjugates, which are
type-(iii) double points of $V_4$ ($d_0^2 + d_1^2 = 0$, two pairs zero). So they
lie on the cone lines over $\mathrm{Sing}(V_4)$, at the places where the branch
divisor $d_0 = \sigma$ crosses them. The local singular scheme is
positive-dimensional there (Singular `vdim` = $-1$): they are special points on a
curve of $cA_1$ singularities.

**Status of "$V_5$ is of general type".**

- The only isolated singularities found are the 128 cone-over-$V_4$ points.
- The rest of $\mathrm{Sing}(V_5)$ consists of curves: the 864 cone lines
  (generic type $cA_1$, canonical) with finitely many special points, and five
  further curve families whose transversal type is pending (B1).
- Any non-canonical locus of dimension $\le 1$ costs $O(m^4)$ conditions, with
  constants bounded by the multiplicities and curve degrees (order $10^4$). That
  is far below $2.9 \cdot 10^8\, m^4$.

So **$V_5$ is of general type** by this count. A formal proof needs explicit
local constants along the curves (B1). **[N19: the count of these points is at
least 77,824, not 128, and $\sigma = 0$ was missing. The margin is about 34.]**

**B1 done** (`v5_transversal.py`, 5 s). For each of the five curve strata, take an
$\mathbb F_p$-point of the stratum ($p = 32003$; for $K = (0,2)$ and
$(0,1,2,3)$, $p = 32009 \equiv 1 \bmod 8$, since those strata live over
$\mathbb Q(i)$ or $\mathbb Q(\zeta_8)$), cut the local model by a random
hyperplane through it, and compute the local singular-scheme length. It is
**1 in all five cases: a transversal ODP, so $cA_1$ along every curve stratum,
which is canonical.** The slow primary-decomposition route
(`v5_local_type.py`) timed out after an hour and is superseded.

**Conclusion.** The non-canonical locus of $V_5$ is finite: the 128
cone-over-$V_4$ points ($a = -3$), plus possibly finitely many special points
where strata meet (e.g. the $K = (0)$ points on the cone lines). Write $\pi$ for
a resolution with exceptional divisors $E_i$ over those points. Then
$K_{\tilde V} = \pi^*\mathcal O(18) + \sum a_i E_i$ with bounded $a_i$.
Vanishing to order about $\lvert a_i\rvert m$ at finitely many points of small
multiplicity costs $O(m^4)$ conditions with constants of order $10^4$, against
$\mathrm{vol}(\mathcal O(18)) \cdot m^4 / 4! \approx 2.9 \cdot 10^8\, m^4$. So
**$V_5$ is of general type**, modulo an explicit discrepancy bound at the
finitely many special points. Bombieri–Lang therefore predicts $V_5(\mathbb Q)$
is not Zariski dense.

### N13. G5: a proven genus-2 descent, and seven rank-1 quotients (`g2descent.gp`, `g2_rank1_certify.gp`)

**The descent.** The genus-2 quotients
$C_S: y^2 = \prod_{e \in S}(t + e)$ of a bottom-lift curve have all six
Weierstrass points rational. In that case 2-descent is linear algebra over
$\mathbb F_2$:

- **odd model:** $W^2 = \prod (X - a_i)$;
- **descent map:** the $x - T$ map into $(\mathbb Q^*/\mathbb Q^{*2})^5$;
- **local images:** sampled until they reach their known sizes, $2^4$ at odd
  $p$, $2^6$ at 2 and $2^2$ over $\mathbb R$;
- **Selmer group:** the kernel of the stacked local conditions, giving
  $\mathrm{rank} \le \dim \mathrm{Sel} - 4$.

Validation (`g2descent_test_all.gp`). There are 56 split Jacobians
$J \sim E_1 \times E_2$ (symmetric $S$) with ranks known from PARI. The results:
**0 violations; bound sharp in 52 of 56.**

**Survey.** For 400 random non-symmetric $S$ over Symes' $E_4$, the rank bounds
are:

| bound | 1 | 2 | 3 | 4 | 5 |
| --- | --- | --- | --- | --- | --- |
| curves | 77 | 158 | 124 | 34 | 7 |

Two runs had a short 2-adic image, which makes the bound unproven, so they are
excluded.

**Certified** (`g2_rank1_certify.gp`). For 7 subsets:

- the descent bound is 1, with every local image full;
- $J(\mathbb Q)_{\rm tors} = J[2]$ (gcd of $\#J(\mathbb F_p)$, $p < 400$ good, is 16);
- the trivial-lift point's image is outside the 2-torsion span, so it has
  infinite order.

**So $\mathrm{rank}\,J(C_S)(\mathbb Q) = 1 < g = 2$ for**
$S = [1,5,6,11,12,14], [6,7,8,9,14,15], [1,9,11,12,14,15], [2,3,8,9,11,13],
[3,5,8,12,13,16], [5,8,9,11,14,15], [3,6,7,8,9,13]$
(indices into Symes' 16 signed roots). Chabauty–Coleman applies to each of these
curves.

### N14. G5: Chabauty script drafted, untested (`chabauty_g2.py`)

This is rank-1, genus-2 Chabauty in the integral odd model, at a good prime
$p \ge 5$:

1. **Annihilating differential.** Compute $\omega$ from one Coleman integral
   $\int_{\iota P_0}^{P_0}$, which is $2\log(P_0 - \infty)$.
2. **Discs of known points** (Weierstrass points, $\infty$, $P_0$,
   $\iota P_0$): tiny integrals only. Their global integrals vanish, being
   2-torsion or annihilated.
3. **Empty discs:** one Coleman integral from $P_0$ to a lift.
4. **Zero counts:** Strassmann bounds per disc.

Certificate iff every disc's zero bound equals its number of known points. The
8 known points are all degenerate lifts, so a certificate means **no $E_5$
above Symes' $E_4$, for any $\mu$**. It needs Sage (passagemath, being
installed).

### N15. G5: the first bottom-lift certificate: no $E_5$ above Symes' $E_4$ (`chabauty_g2.py`, `mw_sieve.py`)

**Smoke test** (`chabauty_smoke.py`, passagemath 10.8.13 in `~/sage_venv`). For
every rational Weierstrass point $W$,
$\int_W^{P_0}\omega = \tfrac12\int_{\iota P_0}^{P_0}\omega$ holds to
$O(17^{15})$. So Sage's Coleman integrals and our odd models agree.

**Chabauty.** At $p_0$ = 17 (or 13), with the known points in distinct discs,
every disc has Strassmann bound exactly 1. So each known disc holds only its
known point, and each empty disc holds exactly one $p$-adic zero. That second
part is forced: integrals of normalized differentials lie in $p\mathbb Z_p$ when
$p \nmid \#J(\mathbb F_p)$, and the integral maps a disc onto $p\mathbb Z_p$.
Chabauty alone therefore never excludes an empty disc. One bug was fixed along
the way: Laurent-series indexing at $\infty$.

**Mordell–Weil sieve.** Here $J(\mathbb Q) = \langle G\rangle \oplus J[2]$ and
$D_0 = [P_0 - \infty] = mG + T_0$, with $m$ odd by the descent.

- **Saturation:** for each odd $\ell$ dividing an order used, $D_0$ is not
  $\ell$-divisible mod $J[2]$. The test: $(\#J(\mathbb F_q)/\ell)(D_0 - T) \ne 0$
  for all 16 $T$, at some prime $q$. Then the image of $J(\mathbb Q)$ in
  $J(\mathbb F_p)$ is $\langle D_0\rangle + J[2]$.
- **The sieve:** an unknown rational point has class $cD_0 + T$, with one integer
  $c$ and one $T$ at every prime. It must land in an empty disc mod $p_0$, and on
  the curve mod the other primes. The search is incremental and greedy: each step
  adds the prime with the smallest expected growth.

**Results** for the 7 rank-1 quotients of Symes' (root-form) lift curve:

| $S$ | $p_0$ | sieve primes | survivors | saturation | |
| --- | --- | --- | --- | --- | --- |
| [1,5,6,11,12,14] | 17 | 67, 23, 103 | 28 → 8 → 2 → **0** | 3, 7, 23 | **certified** |
| [6,7,8,9,14,15] | 17 | — | **0** already at $p_0$ | none needed (ord 8) | **certified** |
| [2,3,8,9,11,13] | 13 | (several) | 4 → **0** | 3, 17 | **certified** |
| [3,5,8,12,13,16] | 17 | (several) | 20 → **0** | 19 | **certified** |
| [5,8,9,11,14,15] | 13 | (several) | 12 → **0** | 3, 11 | **certified** |
| [3,6,7,8,9,13] | 17 | (several) | 20 → **0** | 19 | **certified** |
| [1,9,11,12,14,15] | 17 | pool to 220: …, 181, 193, 101 | 20 → 8 → **0** | 3, 5, 11, 19 | **certified** (default pool: 6 survivors after 12 steps) |

> **Theorem (computational).** The genus-2 curves $C_S$ for all seven certified
> $S$ have exactly 8 rational points, all degenerate lifts ($\mu = 0$, or a zero
> root). Hence **no $E_5$ has Symes' $E_4$ (roots 11, 77, 101, 131, 343, 353,
> 359, 367) as its derived $E_4$, for any scaling $\mu$.** Seven independent
> curves give the same conclusion. The twin (point) form is certified too
> (N17).

### N16. Soundness of the certificate (`omission_test.py`) and the calibration lesson

**Calibration with a real lift is impossible in genus 2.** Every 6-subset
quotient of Symes' derived $E_3$, which genuinely lifts at $t_0 = 67405/12$, has
rank bound $\ge 3$ (`calib_ranks.gp`), and $t_0$ lies on all of them. A real lift
supplies extra independent points, so the rank exceeds the genus. The flip side:
Symes' $E_4$ has many rank-1 quotients (19%), which fits it having no lifts.

**Injection by omission instead.** On certified curve [1,5,6,11,12,14], hide one
genuine rational point (treat its disc mod 17 as empty) and rerun with the same
primes. In all three trials (hiding $\iota P_0$, $W_1$ or $W_3$) the hidden
point's true class $(c, T)$ **survives** and the curve is **not** certified.
Final survivor counts: 1, 2, 1. So the sieve does not eliminate genuine rational
points, and the zero-survivor certificates are meaningful.

### N17. A5: certificates over the catalog (`lift_pipeline.py`, `certify_catalog.py`, running)

The generic pipeline: for any seed, random non-symmetric 6-subsets, then rank-1
certification (`g2_rank1_generic.gp`), then Chabauty, then the MW sieve with
saturation. At about 20–80 s per seed, the results so far are:

| $L_1$ | root form | point form |
| --- | --- | --- |
| 67405 (Symes) | certified, via [2,8,11,12,13,15] (0-indexed), $p_0 = 17$ | certified, $p_0 = 17$ |
| 5915065 | certified, $p_0 = 23$ | certified, $p_0 = 23$ |
| 600445 | certified, $p_0 = 19$ | certified, $p_0 = 19$ |

The 67405 root-form certificate here uses a different quotient from N15's, so it
is an independent re-proof. **Each certified entry: no $E_5$ has that $E_4$ (in
that twin form) as derived $E_4$, for any $\mu$.**

**Completed** (`certify_parallel.py`, 6 workers, about 6 min for the last 27
seeds; logs `certify_catalog.log`, `certify_families.log`,
`certify_parallel.log`). **All 42 seed forms are certified, with no undecided
seed.**

| seeds | forms | certified | Chabauty primes | largest root |
| --- | --- | --- | --- | --- |
| all 15 catalog $E_4$ ($L_1 \le 10^{12}$) | root and point | 30 of 30 | 13–47 | $1.2 \cdot 10^6$ |
| family 1, $k = 2, 3, 4$ (twins: $k = -3, -4, -5$) | root and twin | 6 of 6 | 19–31 | $2.0 \cdot 10^{27}$ |
| family 2, $k = 8, 9, 10$ | root and twin | 6 of 6 | 13, 17 | $3.3 \cdot 10^9$ |

Per seed, the first random non-symmetric rank-1 quotient tried was certified
every time; no seed needed a second quotient.

> **Theorem (computational).** Let $E$ be any of these: an $E_4$ with
> $L_1 \le 10^{12}$; a member of Bremner's family 1 with $k \in \{-5, \dots, 4\}$;
> or a member of family 2 with $5 \le \lvert k\rvert \le 10$ (in either 45° twin
> form). Then **no $E_5$ has $E$, up to scaling, as its derived $E_4$.**
>
> Equivalently: the primitive derived $E_4$ of any $E_5$ lies outside the catalog's
> twin classes and is not one of those family members. This is complementary to
> Proposition 6.1 (no $E_5$ with $L_1 \le 10^{12}$): together, an $E_5$ needs
> $L_1(E_5) > 10^{12}$ **and** $L_2(E_5) > 5 \cdot 10^{11} c^2$, where $c$ is the
> content of the derived sequence.

**Correction (audit, N19).** The first version of this entry said
$L_2(E_5) > 10^{12} c^2$, which overstates the bound by up to a factor of 2.
The catalog's $L_1$ is the paper's $L_1$ of the *root* form. The point form has
odd pair sums and paper $L_1 = L_1^{\rm cat}/2$. So a primitive derived $E_4$
with odd pair sums is covered only up to paper $L_1' \le 5 \cdot 10^{11}$. With
even pair sums it is covered up to $10^{12}$. Since $L_2(E_5) = c^2 L_1'$, the
safe statement is $5 \cdot 10^{11} c^2$. (Family 1 at $k = 0, -1$ is degenerate
and not a seed.)

### N18. C1/D2: the slice and $\varphi$ at larger height (`slice_search.cpp`, `slice_search_gpu.cu`, `phi_scan400.gp`)

- **C++ search** (`slice_search.cpp`, exact `__int128`, 16 threads): height 500 in
  1.1 s, matching the Python search exactly.
- **CUDA residue-table sieve** (`slice_search_gpu.cu`, RTX 4060). For each $z$,
  16 primes $\le 61$ give $p \times p$ bitmasks of residues where $q_1, q_2$ are
  both squares or 0; the GPU filters and the host does the exact test on
  survivors. Hit set identical to the CPU at $H = 500$ and $H = 2000$. Height 2000
  takes 5.8 s, with 82,884 survivors and 78 exact hits.
  - **Build note:** glibc 2.43 on this machine declares `cospi`/`sinpi` as
    `noexcept`, which clashes with CUDA 12.9's headers. Even an empty `.cu` fails.
    `-U_GNU_SOURCE` fixes it. The same issue will hit `code/native`'s CUDA
    preset when it is rebuilt.
- **Heights 2000, 5000 and 10,000** (GPU: 5.8 s, about 90 s, about 12 min;
  12.4 M survivors and 92 exact hits at 10,000): the same 8 slice $E_4$ as at
  1200. Those are the families' members up to family-2 $k = \pm 9$, and the
  exemplar. **Nothing new through height $10^4$.** Family-2 $k = 10$ lies above
  that height.
  - **C1 gate reached:** the exemplar is treated as a sporadic point of $S_B$.
    A third rank-$\ge 1$ curve through it with moderate height growth would
    almost certainly have shown companions by now.
- **D2.** The $\varphi$-partner scan over $\lvert t\rvert, u \le 400$ (about
  240,000 base points, 209 s) finds **0** rational partners with a different
  $E_4$, and 524 special points, all on the $u = \pm 3t \pm 2s$ lines.

## 2026-10-23

### N19. Audit of N1–N18 (`audit_certs.py`, `v5_audit_points.py`; report in `../ag_research_report.md`)

Every PROVEN label was re-derived from its script, and every script-backed claim
re-run or cross-checked. Five findings changed something.

**(a) MW-sieve soundness gap (N15, N17): fixed.**

- **The gap.** The sieve needs the index $m$ of $\langle D_0\rangle + J[2]$ in
  $J(\mathbb Q)$ to be prime to $\#J(\mathbb F_q)$ at every sieve prime. The code
  proved $\ell \nmid m$ only for the odd $\ell \mid \mathrm{ord}(D_0 \bmod q)$, not
  for every odd $\ell \mid \#J(\mathbb F_q)$. N15's docstring states the right
  condition.
- **The fix** (`lift_pipeline.mw_sieve`, on by default). At each $q$, multiply
  every class by $k_q$, the product of the odd $\ell$-parts of
  $\#J(\mathbb F_q)$ with $\ell \nmid \mathrm{ord}(D_0)$. Sound for any $m$:
  with $\Lambda = \mathrm{lcm}_q(\#J(\mathbb F_q)/k_q)$ we have
  $\gcd(m, \Lambda) = 1$, so $k_q[P - \infty] = k_q(cD_0 + T)$ with one $c$ and
  one $T$ for all $q$.
- **The re-run** over all 49 certificates (the 42 seed forms and N15's 7):
  - E4 checks, rank 1 and Chabauty (same $p_0$): 49/49.
  - 16 old certificates used a prime with an untested $\ell$-part. In 9, that
    $\ell$ (3, 5, 7, 19 or 47) was never proven saturated, so the old proof was
    incomplete.
  - The fixed sieve re-certifies 46 of 49 quotients with the default pool.
  - An independent second quotient (seed 2, fixed pipeline) certifies **42/42**
    seed forms.
  - `hyperellratpoints` to height 3000 finds no extra points anywhere.

  - The other 3 quotients (cat12781559905 point form, F1k4 twin form, and N15's
    [1,9,11,12,14,15]) re-certify with the pool widened to 220
    (`audit_retry_wide.py`, `retry3_wide_pool.log`).

  So **N17's theorem stands with the fixed pipeline**: 49/49 quotients, and two
  independent certificates for every one of the 42 seed forms.

**(b) N17's $L_2$ bound** was too strong by up to a factor of 2 (point forms have
paper $L_1 = L_1^{\rm cat}/2$). Corrected in N17. Prop. 6.1 is unaffected, and
the 30 catalog forms have pairwise distinct $(J_1, J_2)$, so no collision exists
across twin forms either.

**(c) $\mathrm{Sing}(V_5)$ is larger than N12/N12a recorded**
(`v5_audit_points.py`):

- **Cone-over-$V_4$ points, $a = -3$, Jacobian rank 7 (corank 4):**
  - the 128 worst points are one sign pattern of $2^8$: **32,768** in all;
  - over the vertex of the base cone (all roots $\pm 1$) there are **32,768**
    more;
  - at $\sigma = 0$ (skipped by N12), points with half the tree zero: **at least
    12,288** more.
- **Other loci at $\sigma = 0$:**
  - curves where two pairs are entirely zero (rank 10, a rank-4 quadric, so
    $cA_1$);
  - singular over the nodes of $V_4$;
  - smooth otherwise (rank 11).

The plurigenus count becomes $\ge 77{,}824 \cdot 108\, m^4 \approx 8.4 \cdot
10^6\, m^4$ against $2.9 \cdot 10^8\, m^4$: a margin of about 34, not
$2 \cdot 10^4$. General type survives the count, but B2 now needs real
discrepancies at the remaining special points.

**(d) G4's scope.** Family 2's 45° twins are not on the slice ($k = 5$–10
checked), so N8d says nothing about them. The inputs check out: every family
member lies on its own $(s:t:u)$-curve, 4 or 2 of its 8 labelings, and never on
the other.

**(e) B1 robustness.** Six random points per curve stratum: 30/30 transversal
ODPs.

### Blocked / waiting

- **G5 final step:** Sage install (owner). The descent no longer needs Magma.
- **G3(a) (Picard number of the K3s)** needs $\#X(\mathbb F_{p^n})$ for $n$ up to
  about 10, which is too slow in Python. It needs a small C++ point counter.
- **G1.2 local types:** `v5_local_type.py` is running (primary decomposition,
  slow on the first stratum).
