# Algebraic geometry for the $E_5$ problem: a plan of attack and a guided tour

*A plan for attacking the problems of `paper.md` with algebraic geometry and
symbolic computation, written for a reader new to the field. Concepts are
introduced only when the problem needs them. Exercises throughout, answers in
§11, reading list in §10. Every number in a worked example was computed by the
scripts in `drafts/ag_pilots/` (§12), not copied.*

Status: 2026-10-21. Companion to `paper.md` (the objects), `enumeration_primer.md`
(the searches S-A/S-E) and `local_model_plan.md` (the local model S-D).

---

## The bottom line first

Four pilot computations, all run today, already change what the paper can say.

1. **No congruence argument can ever rule out $E_5$.** Every depth-$d$ tower
   variety has nondegenerate $p$-adic points for every prime $p$, by a one-line
   construction from Corollary 3.2 (§4). So the "nonexistence branch" of §6.6(c)
   cannot fire. The valuation floor $V_p(k,d)$ of Proposition 4.3 is bounded in $k$
   for every $p$ and $d$. The explicit towers also **settle three depth-5 entries
   of the §4 table**: at $p = 5, 11, 13$ the floor is exactly $4, 4, 2$. For the
   other primes they give a finish line: $V_3 \le 14$, so $k = 15$ decides $p = 3$.
2. **$V_4$ is of general type, and this is now proved.** Its singular locus is
   exactly 864 ordinary double points, all of them degenerate (§2.3). Together
   with the canonical class $\mathcal O(2)$ from adjunction this gives general
   type. §8.1 claimed general type from adjunction alone, which is not enough on
   its own.
3. **Bremner's two families do not pair into an $E_5$.** This is strong evidence,
   not a proof (§6). I counted mod $p$ at six primes. Within each family there is
   never a pair with the same $(L_1:L_2:L_3)$ and a different $L_4$. Across the two
   families such a pair appears once at $p = 10007$ and never at the other five
   primes. A rational coincidence would show up at all but finitely many primes.
   The proof is a defined elimination computation.
4. **Above a fixed $E_4$ there are only finitely many $E_5$, but the cheap
   certificate is closed.** The $E_5$s above a fixed $E_4$ are the rational points
   of one curve of genus $98{,}305$, so there are finitely many for every scaling
   $\mu$ at once (Faltings). But every elliptic quotient I tested has rank $\ge 1$:
   all 132 of them for Symes' $E_4$. This is structural, not bad luck (§7). A
   certificate therefore needs genus-2 Chabauty, which means Magma or Sage; PARI
   cannot do it.

**The strategic picture.** $E_5$ asks whether a 4-fold of general type has a
nondegenerate rational point. That 4-fold has points over $\mathbb R$ and over
every $\mathbb Q_p$, so local methods are exhausted as a route to nonexistence,
and whatever remains is global. The Bombieri–Lang conjecture says the rational
points of $V_4$ and $V_5$ lie on special subvarieties, such as genus-0 and genus-1
curves like Bremner's families. The plan follows that picture:

- G3: find the special subvarieties of $V_4$.
- G4: decide whether any two of them combine into an $E_5$.
- G5: certify the fibers above individual $E_4$s.

---

## 0. The map

| Line | Question | Main tools | Status today | Deliverable |
| --- | --- | --- | --- | --- |
| **G1** Kummer tower | what kind of varieties are $V_n$? | complete intersections, singularities, adjunction | $V_4$ done; $V_5$ open | §8.1 as a theorem |
| **G2** local–global | can congruences decide $E_5$? | Hensel, $p$-adic points | **done: no** | §6.6(c) corrected; §4 table settled at 3 primes |
| **G3** special subvarieties | where do $V_4$'s rational points live? | K3 surfaces, totally tangent curves, symmetry | Bremner's slice identified | new $E_4$ families, or a reason there are none |
| **G4** top-down | can two family members form an $E_5$? | elimination, resultants, heights | mod-$p$ evidence: no | §6.7 made rigorous |
| **G5** bottom-up | which $E_4$ carry an $E_5$ above them? | Faltings, Jacobians, Chabauty | elliptic route closed | per-seed certificates, all $\mu$ |
| **G6** exploratory | torus/tropical view of $W$; Galois view | Newton polytopes, iterated Galois groups | ideas | time-boxed remarks |

§1 and §3 are background. §2 and §4–§8 each cover one line. §9 is the schedule.

---

## 1. A dictionary: from the paper to geometry

### 1.1 Varieties and their points

A **variety** is the solution set of a system of polynomial equations. Our
equations have integer coefficients, so one system can be solved in many rings:
$\mathbb C, \mathbb R, \mathbb Q, \mathbb Z, \mathbb Q_p, \mathbb F_p,
\mathbb Z$ . Geometry studies the complex solutions: their dimension, degree,
singularities and shape. Arithmetic studies the rational or integral solutions
inside them. The organizing principle of Diophantine geometry is that **the
geometry over $\mathbb C$ constrains the arithmetic over $\mathbb Q$**. For
curves, the genus decides everything (§3). In higher dimension, the Kodaira
dimension conjecturally does (§2.4).

**Projective space.** All our equations are homogeneous in the roots. Scaling
every root by $t$ maps solutions to solutions; this is the paper's
$(L_1, L_2, \dots) \mapsto (t^2L_1, t^4L_2, \dots)$. So the natural home is
projective space $\mathbb P^{N-1}$, the set of nonzero $N$-tuples up to scaling.
A rational point of $\mathbb P^{N-1}$ is exactly a primitive integer solution up
to sign, so "rational points" and "primitive integer solutions" are the same
thing here.

### 1.2 The tower varieties

$V_n \subset \mathbb P^{2^{n-1}-1}$ is the set of root vectors
$(r_0, \dots, r_{2^{n-1}-1})$ satisfying the tree conditions of §2:

- $V_3$ is the quadric surface $r_0^2 + r_1^2 = r_2^2 + r_3^2$.
- $V_4$ is Bremner's threefold in $\mathbb P^7$: three quadrics and one quartic.
- $V_5$ is a 4-fold in $\mathbb P^{15}$ cut out by 11 equations.

An $E_n$ solution is a rational point of $V_n$ off the **degenerate locus**
(Proposition 5.4: a root equal to zero, or two roots equal up to sign).

**Dimension.** $V_n$ has dimension $n - 1$ (Exercise 1.2). That is the dimension
of the parameter space $(L_1 : \dots : L_n)$, which is a weighted projective space,
and this is no accident. The map from roots to parameters is finite, and $V_n$ is
the **splitting variety** of the generic $E_n$. Over $\mathbb C$, every $E_n$
splits, and $V_n \to \{(L_1:\dots:L_n)\}$ is a Galois cover whose group is the
iterated wreath product $W_n = C_2 \wr \cdots \wr C_2$, of order $2^{2^n - 1}$.
Over $\mathbb Q$, a parameter vector is an $E_n$ solution exactly when the Galois
group of the specialized polynomial is **trivial**. Hilbert's irreducibility
theorem says the typical specialization has the full group $W_n$, which is the
first sign that solutions are rare.

### 1.3 One variety, many rings

| Ring | What its points mean | Where in the paper |
| --- | --- | --- |
| $\mathbb C$ | the shape of $V_n$ | §8.1 |
| $\mathbb R$ | real towers (they exist) | — |
| $\mathbb F_p$, $\mathbb Z/p^k$, $\mathbb Z_p$ | local solubility, valuation floors | Prop. 4.3, §6.6 |
| $\mathbb Q$ (= primitive $\mathbb Z$) | $E_n$ solutions | everywhere |
| $\mathbb Q(i)$, $\mathbb Z[i]$ | the multiplicative parameterization | §5 |

**Exercise 1.1.** Rewrite $V_3$ as $(r_0 - r_2)(r_0 + r_2) = (r_3 - r_1)(r_3 + r_1)$
and derive a rational parameterization by two pairs $(\alpha:\beta)$, $(\gamma:\delta)$.
This is Dilcher's $E_3$ theorem and the Brahmagupta–Fibonacci identity of §5,
seen as the statement $V_3 \cong \mathbb P^1 \times \mathbb P^1$.

**Exercise 1.2.** Count the tree conditions of $V_n$ level by level and show the
dimension is $n - 1$, assuming the conditions are independent (§2 proves it).

**Exercise 1.3.** Show $|W_n| = 2^{2^n - 1}$ from $|W_n| = 2\,|W_{n-1}|^2$.

---

## 2. Line G1 — the Kummer tower

### 2.1 Square-root covers

A **double cover** is what you get by adjoining one square root: $\{w^2 = f\}$
lies over the space where $f$ lives, two points over each point with $f \ne 0$
and one point over $f = 0$. The locus $f = 0$ is the **branch divisor**. Adjoin
square roots of $m$ functions at once and you get a **Kummer cover** with group
$(\mathbb Z/2)^m$. It is smooth wherever the branch divisors are smooth and cross
transversally ("normal crossings"). It is singular where they do not, and there
the singularities must be examined.

### 2.2 Proposition 3.1 is a geometric statement

Pass to squared coordinates $s_j = r_j^2$. The $L_1$ conditions become linear:
$s_{2k} + s_{2k+1} = \sigma$ for a common $\sigma = 2L_1$. Put
$d_k = s_{2k} - s_{2k+1}$ (this is $2u_k$ in §3's notation). Every higher condition
then becomes the corresponding tree condition on $(d_k)$, one level down. That is
exactly Proposition 3.1. In one sentence:

> **$V_n$ is the Kummer cover of the cone over $V_{n-1}$ given by
> $r_{2k} = \sqrt{(\sigma + d_k)/2}$, $r_{2k+1} = \sqrt{(\sigma - d_k)/2}$,
> branched along the $2^{n-1}$ hypersurfaces $\sigma = \pm d_k$.**

Read upward, the sentence is Proposition 3.1 (the self-similarity). Read downward,
it is Corollary 3.2 (the bottom-lift criterion). Three consequences follow.

**(a) $V_n$ is a complete intersection.** Each squaring step doubles degrees.
Level $j$ contributes $2^{n-1-j} - 1$ equations of degree $2^j$, for
$j = 1, \dots, n-1$:

| $n$ | ambient | type | degree | $\dim$ |
| --- | --- | --- | --- | --- |
| 3 | $\mathbb P^3$ | $(2)$ | 2 | 2 |
| 4 | $\mathbb P^7$ | $(2,2,2,4)$ | 32 | 3 |
| 5 | $\mathbb P^{15}$ | $(2^7, 4^3, 8)$ | $2^{16} = 65{,}536$ | 4 |

**(b) The canonical class.** The canonical class $K$ measures how "curved" a
variety is: negative for rational varieties like $\mathbb P^n$, zero for K3 and
Calabi–Yau varieties, positive for general type. For a complete intersection,
**adjunction** gives $K = \mathcal O(\sum \deg - \#\text{variables})$, which here is
$$K_{V_n} = \mathcal O\big((n-4)\,2^{n-1} + 2\big):\qquad
K_{V_3} = \mathcal O(-2),\quad K_{V_4} = \mathcal O(2),\quad K_{V_5} = \mathcal O(18).$$
This settles the §8 TODO "confirm $V_5$'s multidegree and canonical class".

**(c) The catch.** Adjunction computes $K$ on the smooth part. To conclude
"general type" from $K = \mathcal O(2)$ you also need the singularities to be mild
(*canonical*). Otherwise a resolution can have a smaller canonical class, and the
conclusion fails. §8.1 skipped this step. §2.3 supplies it for $V_4$.

### 2.3 Worked example: the singularities of $V_4$ (computed)

`v4_singular.sing` computes the singular locus of $V_4$ exactly over $\mathbb Q$.
The singular scheme is **864 reduced points**. "Reduced" means each point has
Tjurina number 1, which for a threefold means an **ordinary double point** (ODP):
locally $x^2 + y^2 + z^2 + w^2 = 0$, the mildest singularity there is. Their images
in squared coordinates are 45 points, of four kinds:

| kind (in $s$-coordinates) | base points | points of $V_4$ above each | total | what the roots look like |
| --- | --- | --- | --- | --- |
| all $s_j$ equal (the cone's vertex) | 1 | 128 | 128 | all $\lvert r\rvert$ equal |
| one root of every pair is 0 ($d_k = \pm\sigma$) | 16 | 8 | 128 | four roots 0 |
| $\sigma = 0$, two whole pairs 0 | 12 (8 real) | 8 | 96 | $L_1 = 0$ |
| in each half: one pair has a root 0, the other pair has equal $\lvert r\rvert$ | 16 | 32 | 512 | two roots 0 |

All 864 lie in the degenerate locus. (I first counted these by hand and got 352.
I had missed the last row, where branch divisors from *opposite* halves of the
tree meet along a reducible conic. Compute, don't count.)

> **Theorem G1.1.** $V_4$ is a normal, irreducible complete intersection of type
> $(2,2,2,4)$ whose singular locus is exactly 864 ordinary double points, all
> degenerate. Hence its singularities are terminal and Gorenstein,
> $K_{V_4} = \mathcal O(2)$, and **$V_4$ is of general type**.

*Proof sketch.*

1. A complete intersection is Cohen–Macaulay. Here the singular locus is finite,
   so $V_4$ is normal (Serre's criterion).
2. A positive-dimensional complete intersection is connected, and a normal
   connected variety is irreducible.
3. A threefold ODP is terminal: blowing it up has discrepancy 1.
4. So the plurigenera of a resolution are those of $\mathcal O(2)$, which is
   ample.

$\square$

### 2.4 What general type buys, and what is open

The **Bombieri–Lang conjecture** says the rational points of a variety of general
type are *not Zariski dense*: they lie in a proper closed subset. For $V_4$, that
means every $E_4$ lies on finitely many curves and surfaces, plus finitely many
sporadic points. Bremner's two families are two such curves (genus 1). The
conjecture is not proved for $V_4$, but it is the right map of the territory, and
G3 turns it into a search.

**Open (G1.2): the singularities of $V_5$.** The same analysis is harder one level
up. The cone over $V_4$ is singular along the 864 lines over $V_4$'s double
points. On the lines over the 128 "all $\lvert r\rvert$ equal" points there are
points where 8 branch divisors $\sigma = \pm d_k$ pass at once, and the base is
singular there too. For comparison, 8 branch divisors through a *smooth* point of
a 4-fold would give a cone over a complete intersection of four quadrics in
$\mathbb P^7$. That has $K = \mathcal O(0)$, so the singularity would be log
canonical but *not* canonical. $V_5$ is very probably of general type, since
$K = \mathcal O(18)$ is very positive and the bad points are few, but this needs a
real computation. Do it before §8 claims Lang–Vojta finiteness for $V_5$.

**Exercise 2.1.** With $s_{2k} = (\sigma + d_k)/2$ and $s_{2k+1} = (\sigma - d_k)/2$,
show that the $L_{2,b}$ condition $s_0s_1 + s_2s_3 = s_4s_5 + s_6s_7$ becomes
$d_0^2 + d_1^2 = d_2^2 + d_3^2$. This is Proposition 3.1 for $n = 4$.

**Exercise 2.2.** Derive the type, degree and canonical class of $V_5$ in the table.

**Exercise 2.3.** Show that $(A, \dots, H) = (1,0,1,0,1,0,1,0)$ lies on $V_4$ and
that the Jacobian matrix has rank $< 4$ there. Which row of the table is it in?

**Exercise 2.4.** Why are there 32 points of $V_4$ above each base point in the
last row, and 8 above each base point in the second row?

---

## 3. Curves: the genus decides

### 3.1 Genus and Riemann–Hurwitz

A smooth projective curve over $\mathbb C$ is a surface with $g$ handles; $g$ is
its **genus**. The tool for computing it is **Riemann–Hurwitz**. If $X \to \mathbb P^1$
has degree $N$ and ramification index $e_P$ at each point $P$, then
$$2g(X) - 2 = -2N + \sum_P (e_P - 1).$$
For example, $y^2 = f(t)$ with $\deg f = 2g + 2$ is a double cover branched at
$2g + 2$ points and has genus $g$.

### 3.2 The trichotomy

| genus | rational points | in our problem |
| --- | --- | --- |
| 0 | none, or infinitely many and parameterizable | one congruum: $L \pm u = \square$ (a conic) |
| 1 | none, or a finitely generated group $E(\mathbb Q) \cong \mathbb Z^r \oplus T$ | Bremner's families (rank $r = 1$) |
| $\ge 2$ | **finitely many** (Faltings 1983) | the fibers below |

For genus 1, the **Mordell–Weil theorem** says the rational points form a finitely
generated group; $r$ is the **rank** and $T$ the torsion. **Mazur's theorem**
gives $|T| \le 16$. The rank is hard to compute. Points you find bound it from
below, and **2-descent** bounds it from above. When the two bounds meet, the rank
is proven. PARI's `ellrank` returns both bounds.

### 3.3 The bottom-lift fiber is a curve of huge genus

Fix a seed: an $E_{n-1}$ with $m$ positive roots $\rho_k$. The $E_n$s above it
are the pairs $(L, \mu)$ with all of $L \pm \mu\rho_k$ squares (Corollary 3.2).
Write $a_k^2 = L + \mu\rho_k$ and $b_k^2 = L - \mu\rho_k$. Eliminating $L$ and
$\mu$ leaves $2m - 2$ quadrics in $\mathbb P^{2m-1}$, which define a curve $D$.
It is a Kummer cover of the $(L:\mu)$-line branched at the $2m$ points
$L/\mu = \mp\rho_k$. Riemann–Hurwitz gives (Exercise 3.1)
$$g(D) = 1 + 2^{2m-2}(m-2).$$

| seed $\to$ lift | $m$ | $g(D)$ | consequence |
| --- | --- | --- | --- |
| $E_2 \to E_3$ | 2 | 1 | infinitely many lifts possible (and they occur) |
| $E_3 \to E_4$ | 4 | 129 | finitely many $E_4$ above each $E_3$ |
| $E_4 \to E_5$ | 8 | **98,305** | **finitely many $E_5$ above each $E_4$** |

> **Proposition G5.1 (finiteness per seed).** For any $E_4$ with distinct nonzero
> roots, only finitely many $E_5$, up to scaling, have it (up to scaling) as
> their derived $E_4$. In particular, §6.4's search over $\mu$ looks for points
> on one fixed curve of genus 98,305, and it can only ever find finitely many.

The table also explains a fact from the data. A family of $E_4$ must move its
derived $E_3$, because no single $E_3$ supports infinitely many lifts. Faltings'
theorem is *ineffective*: it says "finite" but not where the points are. G5 is
about making it effective.

**The top fiber is finite too.** Fix $(L_1, L_2, L_3)$ in general position. The
values $\sqrt{L_4}$ that make $E_4$ split require $E_3(x) = s$ and $E_3(x') = -s$
to have rational solutions. That puts $(x, x')$ on the plane octic
$E_3(x) + E_3(x') = 0$, which is smooth of genus 21 for general $(L_1,L_2,L_3)$.
So only finitely many $L_4$ occur.

### 3.4 Jacobians split into pieces you can compute

The Jacobian of $D$ is an abelian variety of dimension 98,305. Up to isogeny it
splits as a product of Jacobians of hyperelliptic curves
$$C_S:\ y^2 = \prod_{e \in S}(t + e),\qquad S \subseteq \{\pm\rho_0, \dots, \pm\rho_7\},\ |S| \text{ even},$$
of genus $|S|/2 - 1$ (Exercise 3.2 checks the genus total). Every rational point
of $D$ maps to a rational point of *every* $C_S$. So one $C_S$ whose rational
points can be listed is enough to bound the lifts. §7 follows this up.

**Exercise 3.1.** Derive $g(D)$ from Riemann–Hurwitz. (The group is
$(\mathbb Z/2)^{2m}/\langle -1\rangle$; each branch point has ramification index 2.)

**Exercise 3.2.** Check $\sum_{k=1}^{8} \binom{16}{2k}(k - 1) = 98{,}305$.

**Exercise 3.3.** Show that $E_3(x) + E_3(x') = 0$ is smooth at infinity and,
for general $(L_1, L_2, L_3)$, smooth everywhere. Deduce genus 21.

---

## 4. Line G2 — local–global: why congruences cannot finish the job

### 4.1 $p$-adic numbers in one paragraph

A **$p$-adic integer** is a compatible sequence of residues mod $p, p^2, p^3, \dots$;
they form the ring $\mathbb Z_p$, with fraction field $\mathbb Q_p$ and valuation
$v_p$. Squares are easy to recognize by **Hensel's lemma**. For odd $p$, a unit is
a square in $\mathbb Z_p$ iff it is a square mod $p$. For $p = 2$, a unit is a
square iff it is $\equiv 1 \pmod 8$. An integer solution gives a $\mathbb Z_p$
solution for every $p$ and a real solution. If any of these **local** conditions
fails, there is no rational point: a *local obstruction*. The **Hasse principle**
(local everywhere $\Rightarrow$ global) holds for quadrics (Hasse–Minkowski) but
fails in general, for example for Selmer's cubic $3x^3 + 4y^3 + 5z^3 = 0$.

### 4.2 The theorem

> **Theorem G2.1 (no local obstruction, at any depth).** For every $d \ge 1$ and
> every prime $p$ there is a nondegenerate $E_d$ over $\mathbb Z_p$, meaning $2^d$
> distinct nonzero roots in $\mathbb Z_p$; likewise over $\mathbb R$.

*Proof.* Induct on $d$, starting from $x^2 - 1$. Let $\rho_k$ be the positive
roots of a nondegenerate $E_{d-1}$ over $\mathbb Z_p$. Take $L = 1$, and take
$\mu = p$ for odd $p$ or $\mu = 8$ for $p = 2$.

- Every $1 \pm \mu\rho_k$ is $\equiv 1 \pmod p$ (resp. $\pmod 8$), hence a unit
  square in $\mathbb Z_p$ (Exercise 4.1).
- If $L'_j$ are the seed's parameters, the identity
  $E_d(x;\, L, \mu^2L'_1, \mu^4L'_2, \dots) = \prod_k (x^2 - L - \mu\rho_k)(x^2 - L + \mu\rho_k)$
  holds over any commutative ring. So these square roots are the $2^d$ roots of an
  $E_d$ over $\mathbb Z_p$.
- The $1 \pm \mu\rho_k$ are pairwise distinct because the $\lvert\rho_k\rvert$ are
  distinct and nonzero, and they are units, so the new roots are distinct and
  nonzero.

Over $\mathbb R$, take $L = 1$ and $\mu < 1/\max\rho_k$. $\square$

For $d = 5$ the seed can be any integer $E_4$ from the catalog. What the
construction does *not* give is integer square roots of $1 \pm p\rho_k$. That is
the point of $D$ in §3.3, and it is a global question.

**Consequences.**

1. **§6.6(c) can never happen.** Reducing such a tower mod $p^k$ shows
   $V_p(k,d) \le \mathrm{ht}_p(\text{tower})$ for every $k$. The floor is bounded,
   and "$V_p \to \infty$, hence nonexistence" is impossible at every prime and
   every depth.
2. **Any proof that no $E_5$ exists must be global**: Brauer–Manin, descent on
   fibers (G5), or finiteness plus enumeration. No congruence can do it.
3. **The §4 table can be settled from above.** Write $m_p(d)$ for the minimum of
   $\mathrm{ht}_p$ over nondegenerate $\mathbb Z_p$-towers. The proof of
   Proposition 4.3 works verbatim over $\mathbb Z_p$ and gives
   $V_p(k,d) \le m_p(d)$ for all $k$. Explicit towers bound $m_p(d)$ from above.
   When the two bounds meet, the entry is final: no local argument can improve it.

### 4.3 Worked example: the depth-5 table, settled at three primes (computed)

`local_upper.py` searches over all 24 seeds (the 12 catalog $E_4$, each in its two
45°-twin forms, defined in §12). It tries every $v_p(\mu) \le 6$ and every
$L < 64$, and keeps the tower with the smallest maximum node valuation.
`witness_check.py` re-checks the winners with an independent tree-walk oracle,
which reproduces the paper's own exemplar floors (6 at $p = 3$, 2 at $p = 7$):

| $p$ | §4 lower bound $V_p(k,5)$ ($k$ reached) | explicit tower: $\max_{\text{nodes}} v_p$ | verdict |
| --- | --- | --- | --- |
| 2 | 13° ($k = 14$) | 29 (600445, root form, $v(\mu) = 2$, $L = 5$) | $13 \le m_2 \le 29$ |
| 3 | 9° ($k = 9$) | 14 (67405, point form, $\mu = 3$, $L = 1$) | $9 \le m_3 \le 14$; **$k = 15$ decides** |
| 5 | 4† ($k = 5$) | **4** (166045477, point form, $\mu = 1$, $L = 5$) | **settled: 4** |
| 7 | 3° ($k = 3$) | 10 (67405, point form, $\mu = 7$, $L = 1$) | $3 \le m_7 \le 10$ |
| 11 | 4° ($k = 4$) | **4** (43937864165, point form, $\mu = 1$, $L = 3$) | **settled: 4** (the ° was not "still rising") |
| 13 | 2† ($k = 3$) | **2** (12781559905, point form, $\mu = 1$, $L = 13$) | **settled: 2** |

"Settled" assumes §4's lower bounds are right. They come from the C++ dynamic
program, which agrees exactly with its brute-force oracle on ten $(p, k, d)$
fixtures (§6.6).

**Inside the $p = 13$ witness.** The seed is the point form of
$L_1 = 12781559905$, with roots $(82796, 76983, 111992, 15471, 101256, 50287, 105836, 39753)$.
Take $L = 13$ and $\mu = 1$.

- Every $13 \pm \rho_k$ is a unit mod 13 and a square there, hence a square in
  $\mathbb Z_{13}$. It works because $13 \equiv 1 \pmod 4$ makes $-1$ a square,
  and all eight roots happen to be quadratic residues mod 13.
- The node valuations by level, leaves to top, are $(0, 0, 0, 1, 2)$. So
  $m_{13}(5) \le 2$, which matches the lower bound.

**What the precision effect looks like (`padic_oracle.py`).** The $p = 3$ witness
has node profile $(0, 4, 2, 6, 14)$. Reduced mod $3^k$, its top node is
$\equiv 0$ for every $k \le 14$, so the oracle reports floor $= k$ at
$k = 8, 9, 10$. That is exactly the "°" symptom §4 records at $p = 3$, produced
by a tower we know exists. §6.6's caution was right, and now the run has a known
finish line.

**What it says about §4's top-node hope.** Each witness attains its maximum at the
top node. At $p = 5$ and $p = 13$ the top valuation equals CMSV's exponent (4 and
2), so top-node localization can gain nothing there. At $p = 11$ the witness has
top valuation 4, so the hoped-for $11^4$ is the most it could give, and it is
still possible. $p = 2$ and $p = 3$ remain open. The "nine orders of magnitude"
of the §4 TODO therefore rest on 2, 3 and 11 alone.

### 4.4 A certificate the C++ code can produce itself

> **Lemma G2.2 (Hensel certificate).** Let $f = E_d(x; L)$ with $L \in \mathbb Z_p^d$,
> and suppose $a_1, \dots, a_{2^d} \in \mathbb Z$ satisfy
> $v_p(f(a_i)) > 2\,v_p(f'(a_i))$ for every $i$, and
> $v_p(a_i - a_j) \le \min(v_p f'(a_i), v_p f'(a_j))$ for $i \ne j$. Then $f$ has
> $2^d$ distinct roots in $\mathbb Z_p$, so $L$ is a nondegenerate $\mathbb Z_p$-tower.

This is a finite check on residues. The dynamic program can run it on any tower it
finds and print "settled" with a certificate, instead of leaving a ° mark. The cost
is precision: $k$ must exceed $2\max_i v_p f'(a_i)$, which can be much more than
$\mathrm{ht}_p$. The explicit constructions above avoid that cost.

### 4.5 The global obstruction that remains

A variety can have points over every $\mathbb Q_p$ and $\mathbb R$ and still none
over $\mathbb Q$. The main known mechanism is the **Brauer–Manin obstruction**.
For curves it is conjectured to be the only one. For varieties of general type
there is no such expectation, and computing the Brauer group of a 4-fold of degree
65,536 is out of reach. **Verdict: not planned.** The tractable global tool is
descent and Chabauty on the fibers $D$ (G5).

**Exercise 4.1.** Prove with Hensel that $1 + pu$ is a square in $\mathbb Z_p$ for
odd $p$ and $u \in \mathbb Z_p$, and that $1 + 8u$ is a square in $\mathbb Z_2$.

**Exercise 4.2.** At $p = 13$ the witness uses $\mu = 1$, but at $p = 3$ the search
needed $\mu = 3$. Explain why $L \equiv 0 \pmod p$ can work only when $p \equiv 1 \pmod 4$,
and why $\mu = 1$ fails at $p = 3$ for these seeds.

**Exercise 4.3.** Reconcile the table: §4 found $V_3(k,5) = k$ at $k = 7, 8, 9$,
and a $\mathbb Z_3$-tower with $\mathrm{ht}_3 = 14$ exists. Why is there no
contradiction, and what must the dynamic program return at $k = 15$?

**Exercise 4.4.** Prove Lemma G2.2 from the strong form of Hensel's lemma: if
$v(f(a)) > 2v(f'(a))$, there is a unique root $\alpha$ with $v(\alpha - a) > v(f'(a))$.

---

## 5. Line G3 — where the rational points of $V_4$ live

### 5.1 The model: Lang's exceptional set

If Bombieri–Lang holds, every $E_4$ lies on a finite union $Z$ of curves and
surfaces, plus finitely many sporadic points. The components of $Z$ are
*special*: genus-0 or genus-1 curves, or surfaces that are not of general type.
Bremner's families are two components. G3 hunts for the rest. It has two goals:
new infinite families of $E_4$ (more candidates for G4), and an explanation of the
"sporadic" catalog members.

### 5.2 A better model of $V_4$: one root per pair

Use Bremner's labels: pairs $(A,B), (C,D), (E,F), (G,H)$, halves
$\{A..D\} \mid \{E..H\}$. **Fix one root from each pair, $A, C, E, G$.** The four
conditions defining $V_4$ are then *linear* in the squares $b = B^2$, $d = D^2$,
$f = F^2$, $h = H^2$:
$$A^2 + b = C^2 + d = E^2 + f = G^2 + h,\qquad A^2 b + C^2 d = E^2 f + G^2 h.$$
The unique solution (Exercise 5.1) is
$$b = \frac{A^2(E^2 + G^2 - C^2) + C^4 - E^4 - G^4}{A^2 + C^2 - E^2 - G^2},\qquad
d = b + A^2 - C^2,\quad f = b + A^2 - E^2,\quad h = b + A^2 - G^2.$$

> **$V_4$ is birational to the $(\mathbb Z/2)^4$-cover of $\mathbb P^3_{(A:C:E:G)}$
> obtained by adjoining $\sqrt b, \sqrt d, \sqrt f, \sqrt h$.** In other words, an
> $E_4$ is a point $(A:C:E:G) \in \mathbb P^3(\mathbb Q)$ at which four explicit
> rational functions are all squares.

**Worked example (Symes' $E_4$).** $(A, C, E, G) = (-101, -77, 343, 367)$.

- The numerator is $-29{,}433{,}642{,}672$ and the denominator is $-236{,}208$.
- So $b = 124{,}609 = 353^2$, $d = 128{,}881 = 359^2$, $f = 17{,}161 = 131^2$ and
  $h = 121 = 11^2$.

This model gives searches a three-parameter space with four square tests, which
suits *structured* searches along planes and curves in $\mathbb P^3$. It is a
different shape from S-A's per-norm enumeration.

### 5.3 Worked example: Bremner's slice is built from K3 surfaces (computed)

Bremner's linear condition $A - C = E - G$ is a plane in $\mathbb P^3$, which he
parameterizes as $A = -s+t$, $C = s+t$, $E = -s+u$, $G = s+u$. On it the
denominator is $2(t^2 - u^2)$, which divides the numerator of $b$, so all four
functions become quadrics (checked with sympy):
$$B^2 = 5s^2 + 2st + u^2,\quad D^2 = 5s^2 - 2st + u^2,\quad F^2 = 5s^2 + 2su + t^2,\quad H^2 = 5s^2 - 2su + t^2.$$
(Check with Symes: $(s,t,u) = (12, -89, 355)$ gives $B^2 = 353^2$.)

Bremner solves the first two conditions rationally by $(x:y:z) \in \mathbb P^2$.
That leaves two conditions, $F^2 = q_1(x,y,z)$ and $H^2 = q_1(x,y,-z)$, where
$q_1$ is his sextic (2-5). So the slice $S_B$ is, up to birational equivalence,
a **double-double cover of the plane branched along two sextics**.
`bremner_sextic.sing` finds:

- $q_1$ is irreducible.
- Its only singularities are three cusps ($A_2$) at $(1:0:0)$, $(0:1:0)$,
  $(1:1:0)$ on the line $z = 0$, and two nodes ($A_1$) at $(\mp\sqrt5 : \pm\sqrt5 : 1)$.
  The total Tjurina number is 8.
- $q_1(x,y,0) = x^2y^2(x-y)^2$, so the line $z = 0$ is tangent to $q_1$ at three
  points (a tritangent).

A double plane branched along a sextic with only simple singularities is a
**K3 surface**. So $X_1 : w^2 = q_1$ is a K3 surface, and so is $X_2$
($X_2 \cong X_1$ via $z \mapsto -z$). The slice $S_B = X_1 \times_{\mathbb P^2} X_2$
has $K = \pi^*(3H)$: general type, provided $q_1$ and $q_2$ meet mildly (not yet
checked).

Why this matters. K3 surfaces sit between rational and general type: their
rational points can be dense, and their curves are governed by a lattice, the
**Néron–Severi group** $\mathrm{NS}$. A family of $E_4$ in this slice is a plane
curve along which *both* $q_1$ and $q_2$ restrict to squares. Equivalently, it
meets both sextics with even multiplicity everywhere (it is **totally tangent** to
both; Exercise 5.3). Family 2 is such a cubic and family 1 such a quartic. A
totally tangent curve splits in two in $X_1$, giving classes in $\mathrm{NS}(X_1)$.
**So the search for families in Bremner's slice is a lattice problem.** The
geometric Picard number is at least 9: the hyperplane class plus the 8 exceptional
curves of $3A_2 + 2A_1$ (the two nodes are defined over $\mathbb Q(\sqrt5)$). The
tritangent probably adds more.

**The G3 program.**

- **(a)** Compute the Picard number $\rho(X_1)$, bounding it above by van Luijk's
  point-counting method and below by explicit curves.
- **(b)** Find the classes of low degree: conics and cubics totally tangent to $q_1$.
  A direct Gröbner computation is feasible for degree $\le 3$.
- **(c)** Keep the curves totally tangent to $q_2$ as well. Using $z \mapsto -z$,
  that means curves $C$ with both $C$ and its mirror totally tangent to $q_1$.
- **(d)** Compute genus and rank of each survivor. Each one with positive rank is
  a new infinite family. Failing that, prove there are none up to degree $D$.
- **(e)** Look at other slices. Bremner's plane is a tangent plane of the quadric
  $A^2 + C^2 = E^2 + G^2$ along which the denominator divides the numerator of $b$.
  Enumerate all planes with that property, modulo the symmetry group of $V_4$
  (§5.4), and analyze each one the same way.
- **(f)** A cheap census. For the 12 catalog $E_4$ and Bremner's Table 1, test
  every small linear relation $\pm X_i \pm X_j \pm X_k \pm X_l = 0$, modulo
  symmetry. If the known solutions cluster on a few slices, the special locus
  becomes visible directly.

### 5.4 The symmetry group

The signed permutations that respect the tree are $2^8$ sign changes times $2^7$
tree swaps, which gives $2^{14}$ up to the overall sign. Add the 45° twist
$(A,B) \mapsto (A+B, A-B)$, which is Bremner's pairing and the paper's "rotation"
(Exercise 5.4). Use this group to reduce every search in G3 and G4.

**Exercise 5.1.** Derive the formula for $b$.

**Exercise 5.2.** Why is a double plane branched along a smooth sextic a K3 surface
($K = 0$)? Why do simple singularities of the sextic not change that?

**Exercise 5.3.** Show that if $q|_C = F^2|_C$ for a plane curve $C$ that is not a
component of $\{q = 0\}$, then every intersection point of $C$ with $\{q = 0\}$
has even multiplicity.

**Exercise 5.4.** Show that the 45° twist maps $V_4$ to itself.

---

## 6. Line G4 — top-down: twins and family intersections

### 6.1 The forgetful map

Let $\pi: V_4 \to (L_1 : L_2 : L_3)$ forget $L_4$. An $E_5$ is two rational points
in the same fiber with different $L_4$ (the Key Property). The right invariants
are $j_1 = L_2/L_1^2$ and $j_2 = L_3/L_1^4$, which do not change under scaling.
But equal $j$ is *necessary, not sufficient*, over $\mathbb Q$: you also need
$L_1'/L_1$ to be a rational square (Exercise 6.1). If it is not, you get an $E_5$
over $\mathbb Q(\sqrt{L_1'/L_1})$ instead.

### 6.2 Curves meeting curves

Suppose $V_4(\mathbb Q)$ lies on finitely many curves $C_i$. Then $E_5$s come from
pairs $(C_i, C_j)$. Their images $\pi(C_i)$ and $\pi(C_j)$ are curves in a surface,
so they meet in finitely many points unless they coincide. There are two checks
per pair:

1. **Twin components.** A whole curve of coincidences inside one family, or
   between two families. This would be a family of $E_5$s.
2. **The finite intersection, and its rational points.**

### 6.3 Worked example: Bremner's families mod $p$ (computed)

`family2.py` and `families.py` first verify symbolically that Bremner's formulas
satisfy all four equations of $V_4$ modulo his curves. Then they reduce both
families mod $p$, group the nondegenerate points by $(j_1, j_2)$, and count the
groups that contain two different $L_4$. Equal $j$ is weaker than an $E_5$, so
this test over-counts, and zero is the stronger statement.

| $p$ | points (fam. 1, fam. 2) | within-family twins | classes hit by both families | … with different $L_4$ |
| --- | --- | --- | --- | --- |
| 1009 | 936, 950 | 0, 0 | 11 | 0 |
| 2003 | 1976, 2004 | 0, 0 | 10 | 0 |
| 3001 | 2960, 2868 | 0, 0 | 12 | 0 |
| 5003 | 5032, 4902 | 0, 0 | 15 | 0 |
| 7001 | 6808, 6894 | 0, 0 | 10 | 0 |
| 10007 | 9968, 9760 | 0, 0 | 15 | **1** |

Family 2 alone at $p = 100003$ gives 100,525 points, 50,237 classes and 0 twins.

How to read it:

1. **No twin component.** Family 2's map is generically 2:1, which is the known
   symmetry $kQ \leftrightarrow -kQ$, and a twin component would produce about $p$
   collisions.
2. **No accidental collisions within a family either.** This suggests the
   singularities of the image curves sit at degenerate points. That is worth
   proving.
3. **The shared classes with the same $L_4$ are the same $E_4$.** Both families
   lie on Bremner's slice, and two curves on a surface meet.
4. **There is no rational cross-family $E_5$.** A rational one would reduce to a
   collision at every prime outside a finite set (Exercise 6.2). One hit at
   $10007$ and none at five other primes means there is none. This is evidence,
   not proof.

### 6.4 Making it a proof

Compute the scheme
$Z = \{(P, Q) \in C_1 \times C_2 : j(P) = j(Q)\}$ by elimination in Singular,
discard the degenerate components, and confirm that $Z$ is finite. Its rational
points have bounded height. Each family curve has rank 1, so $P = mG_1$ and
$Q = nG_2$, and the height bound becomes a finite range of $(m, n)$ to check. Do
the same for $C \times C$, minus the diagonal and the graph of $P \mapsto -P$.
That is one or two Gröbner computations, possibly heavy: the invariants are
rational functions of degree up to about 200 on the curves. After that, every new
family from G3 goes through the same test.

**Exercise 6.1.** Two $E_4$ have equal $(j_1, j_2)$. Show they form an $E_5$ over
$\mathbb Q$ iff $L_1'/L_1 \in \mathbb Q^{\times 2}$, and over $\mathbb Q(\sqrt{L_1'/L_1})$
otherwise.

**Exercise 6.2.** Why would a rational $E_5$ made of two family members appear as
a collision mod $p$ for all but finitely many primes?

---

## 7. Line G5 — bottom-up: fiber curves and Chabauty

### 7.1 The setting

Fix a seed $E_4$ with roots $\rho_0, \dots, \rho_7$ and put $t = L/\mu$. An $E_5$
above it is a $t$ at which all 16 numbers $t \pm \rho_k$ lie in one square class.
That is a rational point of $D$ (genus 98,305), mapping to every $C_S$ of §3.4.
The cheapest possible certificate would be an elliptic $C_S$ ($|S| = 4$) of rank
0: then $C_S(\mathbb Q)$ is the finite torsion group, which gives finitely many
candidate $t$ to check.

### 7.2 Worked example: every elliptic quotient has rank $\ge 1$ (computed)

For Symes' $E_4$ (`lift_ranks.gp`, `lift_triples_check.gp`):

| family of quotients | curves | ranks found | rank 0 |
| --- | --- | --- | --- |
| $y^2 = \prod_{e\in S}(t+e)$, $\lvert S\rvert = 4$ (sampled) | 20 | 1 (9), 2 (6), 3 (5) | none |
| $E1_{ijk}: Y^2 = (X-\rho_i^2)(X-\rho_j^2)(X-\rho_k^2)$ | 56 | 1 (15), 2 (32), 3 (9) | none |
| $E2_{ijk}: v^2 = u(u-\rho_i^2)(u-\rho_j^2)(u-\rho_k^2)$ | 56 | 1–4 | none |

**The reason is structural.** Two points of infinite order are built into the
problem.

- **The trivial lift $\mu = 0$** ($t = \infty$) gives two rational points
  $\infty^\pm$ on every even-degree $C_S$. Their difference has infinite order:
  checked on all 56 $E2$, where it is $(0, c)$ (Exercise 7.1).
- **The seed's own $L_1$.** For every root, $2L_1 - \rho_k^2$ is the square of
  its partner, so $X = 2L_1$ is a rational point on every $E1_{ijk}$. It has
  infinite order on all 56 (Exercise 7.2).

So the rank-0 shortcut is closed by the structure of the problem. Both points
exist for every seed, so trying another seed is unlikely to reopen it. (Only
Symes' seed was computed.)

### 7.3 What works instead

- **Chabauty–Coleman.** If the Jacobian of a curve $C$ has rank smaller than the
  genus, then $C(\mathbb Q)$ is finite *and computable*, via $p$-adic integration.
  The natural targets are the genus-2 quotients ($|S| = 6$, of which there are
  8,008). For non-symmetric $S$ the Jacobian is typically simple. The trivial
  point forces rank $\ge 1$, so we need rank exactly 1. The elliptic data
  (ranks 1–3 are typical) warn that rank 2 may be common. **Experiment B1**
  measures the rank distribution.
- **If genus 2 fails:** genus-3 quotients ($|S| = 8$, need rank $\le 2$); elliptic
  Chabauty over number fields (Bruin); the Mordell–Weil sieve (Bruin–Stoll), which
  combines several quotients and exploits the known trivial point.
- **Tools.** PARI has none of this. Magma (`Chabauty`, `RankBound`,
  `TwoCoverDescent`, `MordellWeilSieve`) is the standard; Sage covers part of it.
- **Calibration.** This is the "synthetic injection" test of §6.8. Run the same
  pipeline on the derived $E_3$ of Symes' $E_4$, where the fiber has genus 129 and
  the known lift $t = 67405/12$ must come out (12 is the content of the derived
  sequence $(50244, 67284, 57204, 61476)$).

### 7.4 What a certificate buys, honestly

A certificate is the sentence "no $E_5$ lies above seed $X$, for any $\mu$". It
replaces §6.4's screen over $\mu \le 10^6$ with a proof, and it applies to seeds
beyond any bound $B$, such as the first $N$ members of each Bremner family. Its
limit: for catalog seeds it excludes, roughly, the $E_5$ with $L_2/c^2 \le B$,
where $c$ is the content of the derived sequence. That is a *different* region from
Proposition 6.1's $L_1 \le B$, not a larger one.

**Exercise 7.1.** On $y^2 = \prod_{i=1}^4 (t - e_i)$, compute the divisor of
$t - e_i$ and show $2\,[(e_i, 0) - \infty^+] = [\infty^- - \infty^+]$. Deduce that
if $\infty^- - \infty^+$ has infinite order, then so does every
$(e_i, 0) - \infty^+$, and the rank is $\ge 1$.

**Exercise 7.2.** Verify that $X = 2L_1$ lies on $E1_{ijk}$ for every triple.

**Exercise 7.3.** Why are only *even* subsets $S$ usable?

---

## 8. Line G6 — exploratory (time-boxed)

### 8.1 The parameterization is a torus, and $W$ has eight terms

Over $\mathbb Q(i)$ the $L_1$ conditions $n_j\bar n_j = n_1 \bar n_1$ are
*binomial*, so they cut out a toric variety. §5's syndrome parameterization
$X_0, \dots, X_7 \mapsto (n_1, \dots, n_4)$ is its monomial parameterization, and
Proposition 5.1 is surjectivity on integral points. In the variables
$(X_t, \bar X_t)$, the polynomial $W$ is very sparse. Since
$A^2B^2 = -(n_1^2 - \bar n_1^2)^2/16$ and the $\lvert n_j\rvert^4$ terms cancel,
$$W = -\tfrac1{16}\big(n_1^4 + n_2^4 - n_3^4 - n_4^4 + \overline{n_1^4 + n_2^4 - n_3^4 - n_4^4}\big),$$
and each $n_j^4$ is a single monomial. So **$W$ is a sum of 8 monomials in 16
variables**. The size the paper complains about ("far too big to be displayed")
comes entirely from expanding it into real coordinates.

That makes **tropical geometry** concrete. For a sum of terms to vanish
$\mathfrak q$-adically, the minimum valuation must be attained at least twice. At
a split prime with exponents $c_j = v_{\mathfrak q}(n_j)$ and $b = v_p(2L_1)$, the
8 valuations are $\{4c_j,\ 4(b - c_j)\}$. One example (Exercise 8.1): if
$p^2 \parallel 2L_1$, then **$p$ cannot divide exactly three of the four $n_j$**.
The same mechanism gives a per-prime sieve on token distributions, which is a
possible route to the canonization lemma and a filter for the syndrome census.
The $q^4$-solve of `enumeration_primer.md` §7 fits this picture: a pattern with
one single-prime bucket is a degree-4 cyclic (Kummer, again) cover of the other
buckets.

### 8.2 The arboreal view

$V_n$ is the splitting variety of a non-autonomous iterated quadratic. The
literature on **arboreal Galois representations** (Odoni; Jones' survey) studies
when such preimage trees have *large* Galois groups; we need the trivial extreme.
Counting results for thin sets (Cohen, Serre) give rigorous but weak upper bounds
on the number of $E_4$ below a height. That could be a sanity bound for §8's
heuristic. Low priority.

**Exercise 8.1.** Prove the three-of-four statement two ways: tropically (exactly
one of the 8 terms has $\mathfrak q$-valuation 0), and elementarily (mod $p$, from
$W = 0$ and $G^2 + H^2 = 2L_1$).

---

## 9. The plan

### 9.1 Phases

| Phase | Work | Tools | Output / gate |
| --- | --- | --- | --- |
| **0** — done today | pilots §2.3, §4.3, §6.3, §7.2 | Singular, PARI, Python | the four results at the top |
| **1** — days | write Theorems G1.1 and G2.1 into the paper; annotate the §4 table; replace §6.6(c); finish lines for S-D | — | draft edits for review (list in §9.3) |
| **2** — 1–2 weeks | G4 rigorous family intersection; G3(a–c) $\mathrm{NS}(X_1)$ and totally tangent conics/cubics; G3(f) census | Singular, Python | a new family, or "no $E_5$ from Bremner's families" as a proposition |
| **3** — needs Magma or Sage | G5: ranks of genus-2 quotients (B1); Chabauty for Symes, then the catalog; calibration on an $E_3$ seed | Magma | per-seed certificates, or a documented failure mode |
| **4** — weeks | G3(e) other slices modulo symmetry; G1.2 singularities of $V_5$; G6 remarks | Singular | a geometry paragraph for §9 of the paper |

### 9.2 Decision gates

- **A G3 family whose $\pi$-image coincides with another family's.** Stop
  everything: that would be a family of $E_5$ candidates, and G4 decides it at
  once.
- **B1 shows the genus-2 ranks are always $\ge 2$.** Move G5 to the Mordell–Weil
  sieve, or drop it.
- **S-D at $p = 3$, $k = 15$ still returns $V = k$, or at $p = 7$, $k = 11$
  returns more than 10.** That contradicts Theorem G2.1 plus an explicit witness,
  so it is a bug in the dynamic program. These are free regression tests, once
  the orbit rewrite of `local_model_plan.md` §5.1 reaches those precisions.
  Likewise, $V_{11}(5,5)$ must be 4, and that one is within reach today.

### 9.3 Proposed edits to `paper.md` (not applied)

- **§4 table.** Mark the depth-5 entries at $p = 5, 11, 13$ as sharp, citing the
  explicit $\mathbb Z_p$-witnesses. Retire the ° at $p = 11$. Add upper bounds
  $m_2 \le 29$, $m_3 \le 14$, $m_7 \le 10$.
- **§4 TODO (top-node localization).** Note that the possible gain is confined to
  $p = 2, 3, 11$.
- **§6.6(c).** Replace the nonexistence branch with Theorem G2.1: there is no
  local obstruction at any depth. S-D's products are forced divisors and
  densities. Add Lemma G2.2 as the stopping certificate.
- **§8.1.** Replace "hence of general type" with Theorem G1.1. State $V_5$'s type
  and $K = \mathcal O(18)$ from the Kummer-tower count, and flag $V_5$'s
  singularities as open (G1.2).
- **§9.** "Everywhere-local solubility at depth 5" is answered (yes, trivially).
  "Geometry of the solution varieties" becomes the G3 program.
- **§6.7.** Add the mod-$p$ evidence and the proof plan of §6.4 above.

---

## 10. Reading guide

Grouped by line. Within each group, read the first item first.

**Foundations (§1–§2).**
- M. Reid, *Undergraduate Algebraic Geometry*, LMS Student Texts 12, CUP 1988.
  Short and conversational: conics, cubics, what a variety is. Read all of it.
- D. Cox, J. Little, D. O'Shea, *Ideals, Varieties, and Algorithms*, 4th ed.,
  Springer UTM 2015. Ch. 1–3 (Gröbner bases, elimination) and Ch. 8 (projective
  geometry). This is the manual for every Singular computation here.
- J. Harris, *Algebraic Geometry: A First Course*, GTM 133, Springer 1992. Degree,
  dimension and Bézout, done concretely.
- R. Hartshorne, *Algebraic Geometry*, GTM 52. Reference only: complete
  intersections (II.8, III Ex. 5.5), Riemann–Hurwitz (IV.2), adjunction.

**Curves and elliptic curves (§3, §7).**
- J. Silverman, J. Tate, *Rational Points on Elliptic Curves*, 2nd ed., UTM 2015.
  Ch. 1–3: the group law, Mordell's theorem, torsion. The gentlest entry.
- J. Silverman, *The Arithmetic of Elliptic Curves*, 2nd ed., GTM 106, 2009.
  Ch. X: descent and computing ranks.
- J. Cremona, *Algorithms for Modular Elliptic Curves*, 2nd ed., CUP 1997.
  2-descent as an algorithm.

**Diophantine geometry: the big picture (§2.4, §3, §4.5).**
- M. Hindry, J. Silverman, *Diophantine Geometry: An Introduction*, GTM 201,
  Springer 2000. Part A (divisors, canonical class), Part E (Faltings' theorem),
  and Part F (open problems, including Lang's conjectures).
- B. Poonen, *Rational Points on Varieties*, GSM 186, AMS 2017. The early
  chapters (points over various rings) and the chapter on cohomological
  obstructions (Brauer–Manin).
- S. Lang, *Number Theory III: Diophantine Geometry*, Encyclopaedia of Math.
  Sciences 60, Springer 1991. A survey of the conjectures.

**$p$-adic numbers (§4).**
- F. Gouvêa, *$p$-adic Numbers: An Introduction*, 2nd ed., Universitext 1997.
  Ch. 3 is Hensel's lemma in every form used here.
- N. Koblitz, *$p$-adic Numbers, $p$-adic Analysis, and Zeta-Functions*, 2nd ed.,
  GTM 58.
- J.-I. Igusa, *An Introduction to the Theory of Local Zeta Functions*, AMS/IP
  2000. For §6.6(a) and §8: rationality of Poincaré series.

**Chabauty and the effective toolbox (§7).**
- W. McCallum, B. Poonen, "The method of Chabauty and Coleman", in *Explicit
  Methods in Number Theory*, Panoramas et Synthèses 36, SMF 2012, 99–117. Read
  this first.
- M. Stoll, "Rational points on curves", *J. Théor. Nombres Bordeaux* 23 (2011),
  257–277. The whole toolbox in one survey.
- N. Bruin, "Chabauty methods using elliptic curves", *J. reine angew. Math.* 562
  (2003), 27–49.
- N. Bruin, M. Stoll, "The Mordell–Weil sieve: proving non-existence of rational
  points on curves", *LMS J. Comput. Math.* 13 (2010), 272–306.

**Singularities, covers, K3 surfaces (§2.3, §5.3).**
- M. Reid, "Young person's guide to canonical singularities", in *Algebraic
  Geometry, Bowdoin 1985*, Proc. Sympos. Pure Math. 46, AMS 1987, 345–414. Why
  ODPs are harmless and the $V_5$ points might not be.
- R. Pardini, "Abelian covers of algebraic varieties", *J. reine angew. Math.* 417
  (1991), 191–213. The invariants of Kummer covers, such as the $K = \pi^*(3H)$
  of §5.3.
- D. Huybrechts, *Lectures on K3 Surfaces*, CUP 2016. Ch. 1 (examples, including
  double planes) and the chapter on the Picard group.
- R. van Luijk, "K3 surfaces with Picard number one and infinitely many rational
  points", *Algebra & Number Theory* 1 (2007), 1–15. The point-counting method
  for G3(a).

**Computation.**
- W. Decker, C. Lossen, *Computing in Algebraic Geometry: A Quick Start using
  SINGULAR*, Springer 2006.
- G.-M. Greuel, G. Pfister, *A Singular Introduction to Commutative Algebra*,
  2nd ed., Springer 2008.
- The PARI/GP manual: `ellrank`, `ellinit`, `elltors`, `ellorder`,
  `hyperellratpoints`. The Magma handbook: Chabauty, Mordell–Weil sieve,
  Jacobians of hyperelliptic curves.

**Galois theory of iterates (§1.2, §8.2).**
- R. W. K. Odoni, "The Galois theory of iterates and composites of polynomials",
  *Proc. London Math. Soc.* (3) 51 (1985), 385–414.
- R. Jones, "Galois representations from pre-image trees: an arboreal survey",
  *Publ. Math. Besançon* (2013), 107–136.
- J.-P. Serre, *Topics in Galois Theory*, 2nd ed., A K Peters 2008. Ch. 3: Hilbert
  irreducibility and thin sets.

**A four-week path.**
- **Week 1:** Reid, plus CLO Ch. 1–3. Exercises 1.1–2.4.
- **Week 2:** Silverman–Tate Ch. 1–3, plus Gouvêa Ch. 3. Exercises 3.1–4.4.
- **Week 3:** Hindry–Silverman Part A and the Lang part of Part F, plus
  McCallum–Poonen. Exercises 5.1–7.3.
- **Week 4:** Huybrechts Ch. 1, then skim Poonen's Brauer–Manin chapter.
  Exercise 8.1, then G3(a) for real.

---

## 11. Answers

**1.1.** $r_0^2 - r_2^2 = r_3^2 - r_1^2$. Set $r_0 - r_2 = \alpha\gamma$,
$r_0 + r_2 = \beta\delta$, $r_3 - r_1 = \alpha\delta$, $r_3 + r_1 = \beta\gamma$.
Both sides equal $\alpha\beta\gamma\delta$, and every point arises (for
$r_0 \ne r_2$, take $\alpha = r_0 - r_2$, $\gamma = 1$, and solve for $\beta, \delta$).
This is a map $\mathbb P^1 \times \mathbb P^1 \to V_3$ that is bilinear in each
factor. In Gaussian language, $x = gh$ and $y = g\bar h$ do the same thing. A
smooth quadric surface with a rational point and square discriminant is always
$\mathbb P^1 \times \mathbb P^1$.

**1.2.** Level $j$ has $2^{n-1-j}$ groups that must agree, so $2^{n-1-j} - 1$
equations. The total is $\sum_{j=1}^{n-1}(2^{n-1-j} - 1) = (2^{n-1} - 1) - (n-1)$.
The ambient dimension is $2^{n-1} - 1$, so the dimension is $n - 1$.

**1.3.** $|W_1| = 2$, and the exponent satisfies $e_n = 2e_{n-1} + 1$, so
$e_n = 2^n - 1$. Hence $|W_4| = 2^{15}$ and $|W_5| = 2^{31}$.

**2.1.** $s_{2k}s_{2k+1} = (\sigma^2 - d_k^2)/4$, so the condition reads
$\sigma^2/2 - (d_0^2 + d_1^2)/4 = \sigma^2/2 - (d_2^2 + d_3^2)/4$.

**2.2.** Level 1 gives 7 quadrics, level 2 gives 3 quartics, level 3 gives 1 octic.
The degree is $2^7 \cdot 4^3 \cdot 8 = 2^{16}$, and $K = \mathcal O(14 + 12 + 8 - 16) = \mathcal O(18)$.

**2.3.** All pair norms are 1, and both sides of the quartic are 0. The three
quadric gradients are $(2,0,-2,0,0,0,0,0)$, $(2,0,0,0,-2,0,0,0)$ and
$(2,0,0,0,0,0,-2,0)$. The quartic's gradient vanishes because $B = D = F = H = 0$.
So the rank is 3 and the point is singular. It is in the second row: $\sigma = 1$
and $d = (1,1,1,1)$.

**2.4.** In the last row, 2 of the 8 squared coordinates vanish. The other 6 roots
have free signs, modulo the overall sign: $2^6/2 = 32$. In the second row 4 vanish,
giving $2^4/2 = 8$.

**3.1.** The degree is $2^{2m-1}$. There are $2m$ branch points, and over each one
$2^{2m-2}$ points have $e = 2$. So
$2g - 2 = -2 \cdot 2^{2m-1} + 2m \cdot 2^{2m-2} = 2^{2m-1}(m - 2)$.

**3.2.** $\sum_k k\binom{16}{2k} = \tfrac12 \cdot 16 \cdot 2^{14} = 2^{17}$ and
$\sum_{k\ge1}\binom{16}{2k} = 2^{15} - 1$. The difference is
$131{,}072 - 32{,}767 = 98{,}305$.

**3.3.** At infinity the curve is $x^8 + x'^8 = 0$: eight distinct points, so it
is smooth there. An affine singular point needs $E_3'(x) = E_3'(x') = 0$, so $x$
and $x'$ are critical points, and $E_3(x) = -E_3(x')$ means two critical values
are negatives of each other. The critical values are $(L_1^2 - L_2)^2 - L_3$,
$L_2^2 - L_3$ and $-L_3$, and for general parameters no two are negatives. A
smooth plane octic has genus $7 \cdot 6/2 = 21$.

**4.1.** For $f(x) = x^2 - (1 + pu)$ at $x_0 = 1$: $v(f(1)) \ge 1 > 0 = 2v(f'(1))$
when $p$ is odd. For $p = 2$, $v(f'(1)) = 1$, so you need $v(f(1)) \ge 3$, which
is $8 \mid 8u$.

**4.2.** With $L \equiv 0$, you need $\pm\rho_k$ both to be squares mod $p$. That
requires $-1$ to be a square, which means $p \equiv 1 \pmod 4$. At $p = 3$, the
only nonzero square is 1, so $L \pm \rho \equiv 1$ forces $\rho \equiv 0$ for
every root. A primitive seed cannot do that, and the search found no
non-unit alternative with $\mu = 1$.

**4.3.** Mod $3^k$ with $k \le 14$, the witness's top node is $\equiv 0$, so its
capped valuation is $k$. The witness only appears at cap $k$, which looks
"degenerate". Proposition 4.3's floor of 9 is consistent with that. At $k = 15$
the witness has capped profile $\le 14$, so $V_3(15, 5) \le 14$ must hold.

**4.4.** Hensel gives a unique root $\alpha_i$ with
$v(\alpha_i - a_i) > v(f'(a_i))$. Suppose $\alpha_i = \alpha_j$ for some $i \ne j$.
Then $v(a_i - a_j) > \min(v f'(a_i), v f'(a_j))$, which contradicts the
hypothesis. So the $2^d$ roots are distinct, and a polynomial of degree $2^d$ with
$2^d$ distinct roots splits.

**5.1.** Substitute $d, f, h$ (as $b$ plus constants) into the quartic condition.
It is linear in $b$, with coefficient $A^2 + C^2 - E^2 - G^2$.

**5.2.** $K_X = \pi^*(K_{\mathbb P^2} + \tfrac12 \cdot 6H) = \pi^*(-3H + 3H) = 0$,
and $X$ is simply connected. A simple curve singularity gives a rational double
point on the cover. Its minimal resolution is crepant, so $K$ stays $0$.

**5.3.** The divisor of $q|_C$ on (the normalization of) $C$ is the intersection
cycle $C \cdot \{q = 0\}$. It equals $2\,\mathrm{div}(F|_C)$.

**5.4.** $(A + B)^2 + (A - B)^2 = 2(A^2 + B^2)$, so the $L_1$ conditions are
preserved, with $L_1$ doubled. And $(A+B)^2(A-B)^2 = (A^2 - B^2)^2 = (A^2 + B^2)^2 - 4A^2B^2$:
with equal pair norms, the $L_{2,b}$ condition maps to itself.

**6.1.** If $L_1' = \lambda^2 L_1$ with $\lambda \in \mathbb Q$, scale the first
$E_4$'s roots by $\lambda$. Then all three of $L_1, L_2, L_3$ agree, and the Key
Property applies. Otherwise the scaling factor is $\sqrt{L_1'/L_1}$, the scaled
roots are not rational, and the Key Property still holds over
$\mathbb Q(\sqrt{L_1'/L_1})$.

**6.2.** Reduction mod $p$ of the two points, and of the nonzero integers
$L_4 - L_4'$ (numerators) and $\Delta$, is defined and nonzero for all $p$ outside
a finite set.

**7.1.** $\mathrm{div}(t - e_i) = 2(e_i, 0) - \infty^+ - \infty^-$. With origin
$\infty^+$ this says $2P_i = [\infty^-]$. If $[\infty^-]$ has infinite order, so
does each $P_i$, and $r \ge 1$. (Conversely, $r = 0$ would force $\infty^-$ to be
torsion, a special "Pell-type" condition on the $e_i$.)

**7.2.** $2L_1 - \rho_k^2 = \rho_{k'}^2$, where $k'$ is the partner of $k$ in its
pair. So $Y^2 = (\rho_{i'}\rho_{j'}\rho_{k'})^2$.

**7.3.** $t + e = (L + \mu e)/\mu = a^2/\mu$, so all 16 numbers lie in the square
class of the unknown $\mu$. An even product is a square whatever $\mu$ is. An odd
product is $\mu$ times a square, a twist that depends on $\mu$.

**8.1.** *Tropically:* with $c_j = 1$ for three indices and $c_4 \in \{0, 2\}$, the
valuations are $4, 4$ (three times) and $\{0, 8\}$. The minimum 0 is attained once,
so the sum cannot vanish. *Elementarily:* if $p \mid A, \dots, F$, then
$p^4 \mid G^2H^2$, so $p \mid G$ or $p \mid H$. Then $G^2 + H^2 \equiv 0 \pmod p$
forces both, so $p \mid n_4$.

---

## 12. Reproduce

All scripts are in `drafts/ag_pilots/`; run them from that directory. They need
Singular 4.4, PARI/GP 2.17 and Python with sympy and numpy. They read
`catalog_1e11.csv`, the 12-entry S-A catalog to $10^{11}$, regenerated for this
document with `enum_campaign --bound=100000000000
--factor-base=data/factorbase_1e12.bin` (41 s on 16 threads, 10/10 known
solutions recovered). Each catalog row gives four points $z = m + ni$ of norm $L_1$.
They yield two 45°-twin $E_4$:

- the **root form**, with roots $m \pm n$ (integers, $r^2 + s^2 = 2L_1$);
- the **point form**, with roots $m, n$ (half-integral nodes).

| command | what it shows | time |
| --- | --- | --- |
| `Singular -q v4_singular.sing` | $\deg V_4 = 32$; 864 reduced singular points; 45 base points | 1.5 s |
| `Singular -q bremner_sextic.sing` | $q_1$ irreducible; $3A_2 + 2A_1$; total Tjurina 8 | < 1 s |
| `gp -q lift_ranks.gp` | ranks of 20 quartic and 112 triple quotients for Symes' $E_4$ | 1 s |
| `gp -q lift_triples_check.gp` | the undetermined ranks resolved; tallies; trivial point non-torsion on all 56 $E2$ | 4 s |
| `python3 family2.py 1009 10007 100003` | family 2 on $V_4$ symbolically; census | 2 s |
| `python3 families.py 1009 2003 3001 5003 7001 10007` | family 1 symbolically; both families; the §6.3 table | 56 s |
| `python3 padic_oracle.py` | oracle calibrated on the exemplar (6, 2); the $p = 3$ precision effect | < 1 s |
| `python3 local_upper.py` | the §4.3 upper bounds | < 1 s |
| `python3 witness_check.py` | the oracle confirms the $p = 5, 11, 13$ witnesses | < 1 s |

**One cautionary note, in the spirit of §6.8.** The first version of
`local_upper.py` used integer division on the half-integral point-form nodes. It
"beat" the paper's lower bounds at four primes, which is impossible. The
independent oracle caught the error before it reached this document. Keep the
oracle in the loop for every local claim.
