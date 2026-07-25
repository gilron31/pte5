# Split Iterated-Squaring Polynomials and Norm-like Gaussian Integers

Gil Ronen[^1], Roi Garnett[^2]

[^1]: *<gilron31@gmail.com>*
[^2]: *<roigarnett4@gmail.com>*

<!--
STRUCTURE NOTE (restructure in progress):
This paper was restructured from the original research proposal into a
"foundations + experiments" paper, organized BOUNDS-FORWARD:
  1. Introduction
  2. The E_n problem and its PTE correspondence
  3. The tower recursion (the engine)
  4. Divisibility ladder and E_5 lower bounds  [HEADLINE]
  5. Gaussian parameterization and the single-constraint reduction
  6. Algorithms
  7. Computational results — the exhaustive campaign
  8. A density heuristic
  9. Future directions
  Appendix
Inline **TODO:** markers flag gaps that must be filled by the theory/code/
literature workstreams (see TODO.md). Do NOT treat TODO'd claims as proven.
-->

## Abstract

> **TODO:** Rewrite the abstract bounds-forward to match the new structure.
> Lead with: (a) the first explicit lower bounds on any $E_5$ over $\mathbb{Z}$
> (via the tower divisibility ladder), (b) the structure theory that yields
> them, (c) the exhaustive campaign excluding $E_5$ in explicit regions.
> Draft available in `drafts/critical_review.md`. Text below is the original
> proposal abstract, kept temporarily.

We investigate the search for integers $(L_1,L_2,L_3,L_4,L_5)$ such that the polynomial $((((x^2 - L_1)^2 - L_2)^2 - L_3)^2 - L_4)^2 - L_5$ completely factors over the integers. Such a finding will provide an ideal symmetric solution for the Prouhet-Tarry-Escott (PTE) problem with degree $n=15$, extending the state of the art, currently at $n=11$. By linking the problem to the study of norm-like Gaussian integers, we reduced the search problem to finding zeros of a single multivariate polynomial of a special form. In addition, we provide a new enumeration technique for possible solutions. We believe that further results can be achieved following the proposed directions, as well as by running extensive enumerations and examining numerical data.

## 1. Introduction

> **TODO:** Write the introduction. Required beats (see `drafts/critical_review.md`):
>
> - PTE background and state of the art: ideal symmetric solutions known at
>   sizes $\le 10$ and $12$, open at $11$ and $\ge 13$.
> - CMSV [4]: exhaustive symmetric searches at sizes 13–16 (nothing found; size
>   16 only to height 850), and their new ideal solutions over $\mathbb{Z}[i]$
>   and Eisenstein rings at sizes up to 12. Frame this paper as *engaging* CMSV.
> - Gems / BCSS motivation from BMR [1]: an $E_5$ solution is exactly a skew
>   $32$-gem.
> - The recent $E_3$-level density result of Tsai et al. [5] (our §8 starts at
>   $E_4$ because $E_3$ is done).
> - Contributions, stated **bounds-first**: (a) the first explicit lower bounds
>   for any $E_5$ over $\mathbb{Z}$; (b) the tower structure theory producing
>   them; (c) an exhaustive campaign excluding $E_5$ in explicit regions.
> - Identity discipline: this is the $E$-ansatz paper, not a
>   state-of-the-art-at-size-16 paper.
> - Roadmap.

## 2. The $E_n$ Problem and its PTE Correspondence

### The $E_n$ Problem Statement

Given integers $L_1, L_2, \dots, L_n$, we define the associated **$E_n$ Polynomial** recursively as:

$$E_1(x; L) := x^2 - L$$
$$E_n(x; L_1, \dots, L_n) := E_1(E_{n-1}(x; L_1, \dots, L_{n-1}); L_n)$$

The **$E_n$ Problem** is to find integers $(L_1, L_2, \dots, L_n)$ such that the polynomial $E_n(x; L_1, \dots, L_n)$ completely factors over the integers. We refer to such a tuple as an **$E_n$ solution**. We can therefore formulate our primary goal as finding a **non-trivial $E_5$ solution**.

> **TODO (corrected equivalence — this replaces the original overclaim):** The
> original proposal stated that $E_n$ solutions are *equivalent* to ideal
> symmetric PTE solutions of size $2^{n-1}$. The correct statement (per
> `drafts/critical_review.md` §3): a *pair* of $E_n$ solutions sharing
> $L_1,\dots,L_{n-1}$ is equivalent to such a solution (the Key Property below,
> read correctly), but a symmetric ideal solution need not carry the tower
> structure. $E_n$ is therefore a *proper* "recursively symmetric" subclass.
> State this precisely and prove the inclusion is proper via Tarry's 1913
> degree-7 solution $\{0,4,9,23,27,41,46,50\}/\{1,2,11,20,30,39,48,49\}$, whose
> centered positive halves fail the $L_1$ (sum-of-two-squares) condition, so it
> is not an $E_4$. See [1] for the connection between the $E_n$ formulation and
> PTE.

### Necessary and Sufficient Conditions for an $E_n$ Solution

It is clear that if $L_n$ is not a perfect square, then $E_n$ polynomial cannot have integer roots. If $L_n$ **is** a perfect square, we can split $E_n$ polynomial into a product of two $E_{n-1}$ polynomials:

$$
E_n(x; L_1, \dots, L_n) = E_n(x; L_1, \dots, L_{n-1} - \sqrt{L_n})E_n(x; L_1, \dots, L_{n-1} + \sqrt{L_n})
$$

We can continue with this reasoning in a binary tree-like fashion until we recover the full decomposition of $E_n$:

$$
\begin{gathered}
E_n(x; L_1, \dots, L_n) = \prod _{\nu = 0} ^{2^{n-1}-1} (x - r_{\nu})(x + r_{\nu}) \\
r_{\nu} := \sqrt{L_1 \pm \sqrt{L_2 \pm \sqrt{... \pm \sqrt{L_n}}}} \\
\end{gathered}
$$

With $r_0, \dots, r_{2^{n-1}-1}$ being the $2^{n-1}$ positive roots of $E_n$ and the $\pm$ pattern follows the binary expansion of the index $\nu$. The positive roots are canonically sorted, and the constraint $\sqrt{L_k \pm \sqrt{L_{k+1} \pm \sqrt{... \pm \sqrt{L_n}}}}\in \mathbb{Z}$ must hold at every node of the tree. One can rewrite the constraints in terms of the roots. The lowest three constraints are:

$$
\begin{gathered}
r_{2k}^2 + r_{2k+1}^2 = 2L_1 \\
(r_{4k}^2 - r_{4k+1}^2)^2 + (r_{4k+2}^2 - r_{4k+3}^2)^2 = 8L_2 \\
((r_{8k}^2 - r_{8k+1}^2)^2 - (r_{8k+2}^2 - r_{8k+3}^2)^2)^2 + ((r_{8k+4}^2 - r_{8k+5}^2)^2 - (r_{8k+6}^2 - r_{8k+7}^2)^2)^2 = 32L_3
\end{gathered}
$$

We refer to a constraint at the $n$-th level as the **$L_n$ condition**. Note that an $L_{n>1}$ condition possesses certain degrees of freedom arising from the lower conditions in the tree. In [1], an equivalent set of conditions is introduced under the name **Litter Conditions**. Useful augmentations of the $L_2$ condition include:

$$
\begin{gathered}
L_{2,a} := r_0^4 + r_1^4 + r_2^4 + r_3^4 = 4(L_1^2 + L_2)  \\
L_{2,b} := r_0^2r_1^2 + r_2^2r_3^2 = 2(L_1^2 - L_2)  \\
L_{2,c} := (r_0^2 - L_1)^2 + (r_2^2 - L_1)^2 = 2L_2  \\
L_{2,d} := (r_0^2 - r_1^2)^2 - 4r_0^2r_1^2 + (r_2^2 - r_3^2)^2 - 4r_2^2r_3^2 = 8(2L_2 - L_1^2)  \\
\end{gathered}
$$

> **TODO:** Trace the etymology of the "Litter" conditions in [1] and confirm
> our augmentations $L_{2,a\text{-}d}$ against theirs.

### Key Property - Solution Construction

If two $E_n$ solutions share identical $L_1, \dots, L_{n-1}$ coefficients, they can be used to construct an $E_{n+1}$ solution. This suggests a strategy for finding $E_5$ solutions by first finding many $E_4$ solutions, and then find a pair of solutions that share $L_1$, $L_2$ and $L_3$.

> **TODO:** State and prove the construction ($E_{n-1}$-pair $\Rightarrow E_n$)
> as a proposition; this is the "top-split" dual of the bottom-lift criterion
> in §3.

## 3. The Tower Recursion

> **TODO (new section — the structural engine for the bounds in §4):** Develop
> the self-similarity recursion (see `drafts/critical_review.md` §4(i)).
>
> - Define the derived sequence $u_k := (r_{2k}^2 - r_{2k+1}^2)/2$.
> - **Theorem (self-similarity):** the level-$2..n$ conditions on $r$ are
>   exactly the level-$1..(n-1)$ conditions on $u$, up to powers of $2$ in the
>   $L$-normalizations. Hence an $E_n$ solution is precisely an $E_{n-1}$
>   solution in $u$-space, all of whose entries are "congrua" of a common $L_1$
>   ($L_1 \pm u_k$ all perfect squares). Include explicit $2$-power bookkeeping.
> - **Corollary (parity):** with the odd-$L_1$ normalization, $u_k \equiv 0
>   \pmod 4$.
> - **Corollary (bottom-lift criterion):** an $E_5$ exists over a given $E_4$
>   with roots $r'_k$ iff some scaling $\mu$ and some $L$ satisfy $L \pm \mu r'_k
>   = \text{squares}$ for all eight roots — an $8$-way divisor/hash-join. This is
>   a second search axis (lifts from the bottom) dual to the Key Property
>   (splits at the top), and it exposes all meet-in-the-middle cut points.
> - **Novelty check (blocking):** verify against BMR's litter conditions [1] —
>   if they already articulate this level-wise structure, reframe as a
>   congruum-lift reading of [1] rather than "we introduce".

## 4. The Divisibility Ladder and $E_5$ Lower Bounds

> **TODO (new section — the HEADLINE result):** Develop the divisibility ladder
> and the explicit lower bounds (see `drafts/critical_review.md` §4(ii)).
>
> - **Theorem (ladder):** every internal node of the tower is itself a symmetric
>   ideal PTE pair one level down, so the CMSV [4] required divisors of the PTE
>   constant cascade through the nested radicals. For $E_4$: $18 \mid
>   \sqrt{L_3 \pm \sqrt{L_4}}$ and $2^5\cdot3^3\cdot5^2\cdot7^2\cdot11\cdot13 =
>   151{,}351{,}200 \mid \sqrt{L_4}$.
> - **Corollary (first explicit $E_5$ bounds):** using CMSV's $C'_{16}$,
>   any $E_5$ over $\mathbb{Z}$ satisfies (targets to re-derive independently,
>   do NOT trust the review's arithmetic): $\sqrt{L_5} \ge 1.43\times10^{26}$,
>   $L_5 \ge 2\times10^{52}$, largest root $> 3.4\times10^3$, and
>   $L_1 > 5.8\times10^6$ — strictly stronger than CMSV's height-850 exclusion,
>   derived purely from local theory.
> - Verified instances on the catalog exemplar (§7 / Appendix): $v_0,v_1 \equiv
>   0 \pmod{18}$; the smooth factorization of $\sqrt{L_4}$.
> - The forced prime content at every level doubles as a congruence sieve for
>   the search (§6).

## 5. Gaussian Parameterization and the Single-Constraint Reduction

Since the $L_1$ condition for all root pairs implies $r_{2k}^2 + r_{2k+1}^2 = 2L_1$, we can treat the root pairs as Gaussian integers of equal norm. When searching for roots that satisfy higher order $L_n$ conditions, we can confine our search to the complex circle of radius $2L_1$. By fixing $L_1$ and examining its factorization, we can:

1. Determine the number of points with norm $2L_1$ (namely the sum of two squares function, $r_2(n)$ [2]).
2. Efficiently compute the set of all such points. Such a set provides the basis for a search that can yield $E_3$, $E_4$ and higher order solutions.
3. Parameterize sets of norm-like Gaussian integers even without a concrete $L_1$ at hand.

### Factorization Analysis

*The following is the standard derivation of the factorization of an arbitrary rational integer over the Gaussian integers.*

Let $n$ be a norm factoring as:

$$
n = 2^{a_0}p_1^{2a_1} \dots p_r^{2a_r}q_1^{b_1} \dots q_s^{b_s}
$$

- The $q_i$ are odd rational primes such that $q_i \equiv 1 \pmod 4$. They decompose as $q_i = \mathrm{q}_i\overline{\mathrm{q}}_i$.
- The $p_i$ are odd rational primes such that $p_i \equiv 3 \pmod 4$. These remain inert in the Gaussian integers and must appear with even multiplicity. Their existence implies a common factor to all roots. Thus, we generally assume the norm is free of such primes.
- The factor of 2 contributes either a global scaling or a $45^\circ$ rotation, neither of which changes the essential structure of the roots.

Any Gaussian integer $x = a + ib$ with norm $a^2 + b^2 = n$ satisfies (up to units):

$$
x = \prod _{i=1}^s \mathrm{q}_i^{e_i} \overline{\mathrm{q}}_i^{b_i - e_i}
$$

This is characterized by a vector $\mathbf{e} = (e_1, \dots, e_s)$ where $0 \le e_i \le b_i$.

For the sake of simplicity, in the next sections we only treat norms without multiplicities. This allows us to view $\mathbf{e}$ as a binary vector. In that setting, taking the logical NOT of the $e_i$ is equivalent to taking the complex conjugate of $\mathrm{q}_i$.

> **TODO (surjectivity — the load-bearing lemma for every exclusion claim):**
> The multiplicity-free restriction above is exactly the gap that makes any
> "exhaustive to $B$" claim unprovable (see `drafts/critical_review.md` §"failure
> modes", S3): a norm with a square factor $p^2$, $p\equiv1\pmod4$, admits
> *primitive* solutions (one root pair uses the $q^2$ representation, another
> $q\bar q$), and $\sim39\%$ of integers are non-squarefree. Extend the
> parameterization to all admissible norms — exponent vectors $0\le e_i\le b_i$
> (including the real factor $q\bar q = p$), unit twists, and the ramified prime
> $2$ — and prove surjectivity. This is a joint math+code obligation (§6).

### Joint Parameterization of $n$ Norm-like Gaussian Integers

We show it is possible to parameterize any set of $n$ Gaussian integers sharing the same norm using $2^{n}$ variables. We are interested in essentially different solutions, therefore we pay attention to common factors and the order of the roots.

#### Gaussian Integer Pairs and Application to the $E_3$ Problem

Let $x=a+ib, y=c+id$ be norm-like Gaussian integers associated with vectors $\mathbf{e_x}$ and $\mathbf{e_y}$. Let $g = m + in$ be the product of factors where $\mathbf{e_x}, \mathbf{e_y}$ are identical and, $h = p + iq$ where they differ.

$$
\begin{gathered}
x = gh = mp - nq + i(np - mq) \\
y = g\overline{h} = mp + nq + i(np + mq) \\
\end{gathered}
$$

This yields the Brahmagupta–Fibonacci parameterization of $a^2 + b^2 = c^2 + d^2$. Since the $L_1$ condition is the only constraint for the $E_3$ problem, every instance of this parameterization is a valid $E_3$ solution, as correctly pointed out in [3]. The converse also holds. Every $E_3$ solution can be parameterized by  some $g,\ h$.

#### Gaussian Integer Quadruplets and Application to the $E_4$ Problem

To parameterize four norm-like Gaussian integers $n_1, n_2, n_3, n_4$, we identify $2^3 = 8$ distinct "syndromes" (subsets of indices where factors are conjugated relative to $n_1$). Let $X_i = s_i + it_i$ be the product of all factors matching the $i$-th syndrome. Following this observation, we can write down the four Gaussian integers as:

$$
\begin{gathered}
n_1 = A + iB = X_0 X_1 X_2 X_3 X_4 X_5 X_6 X_7 \\
n_2 = C + iD = X_0 X_1 X_2 X_3 \overline{X_4 X_5 X_6 X_7} \\
n_3 = E + iF = X_0 X_1 \overline{X_2 X_3} X_4 X_5 \overline{X_6 X_7}\\
n_4 = G + iH = X_0 \overline{X_1} X_2 \overline{X_3} X_4 \overline{X_5} X_6 \overline{X_7} \\
\end{gathered}
$$

Each of the coefficients $A,B,C,D,E,F,G,H$ can be thought of as a degree 8 polynomial in 16 variables $\bold{X} = \{s_i + it_i\}_{i=0}^{7}$.

#### The General Case

This generalizes to any number of integers. The exponential ($2^{n-1}$) growth of the number of possible syndromes rapidly makes the parameterizations hard to utilize. For e.g., in light of the $E_5$ problem, we can parameterize an octet of norm-like Gaussian integers. The expressions for the $L_2$ and $L_3$ conditions would result in dense polynomials with $2^8 = 256$ variables of degree no less than $512$ and $1024$ respectively.

### Formulation of the $E_4$ Problem Using a Single Constraint

The $L_2$ condition can be written as (Here is the $L_{2,b}$ form, but the result is agnostic to this augmentation):

$$
\begin{gathered}
W(\bold{X}) := A(\bold{X})^2 B(\bold{X})^2 +  C(\bold{X})^2 D(\bold{X})^2 - E(\bold{X})^2 F(\bold{X})^2 - G(\bold{X})^2 H(\bold{X})^2 = 0\\
\end{gathered}
$$

By the construction of the parameterization, the $L_1$ condition:

$$
\begin{gathered}
A(\bold{X})^2 + B(\bold{X})^2 = C(\bold{X})^2 + D(\bold{X})^2 \\
= E(\bold{X})^2 + F(\bold{X})^2 = G(\bold{X})^2 + H(\bold{X})^2\\
\end{gathered}
$$

is automatically satisfied. Therefore the $E_4$ problem now translates to finding the zeros of $W(\bold{X})$, a degree 32 polynomial in 16 variables.

> **TODO:** Hedge the "first single-constraint reduction" claim (per
> `drafts/critical_review.md` §3): symmetric-function elimination also yields
> one resultant equation; the real contribution is a single *polynomial*
> constraint inside an explicit multiplicative parameterization.

Investigating $W(\bold{X})$ is a difficult task due to its size. The actual polynomial expression is far too big to be displayed here, and not even easily manipulated with a computer program. The following observation can be used to simplify the expression: Our quadruplet will not always contain all possible syndromes $X_0,  \dots  , X_7$, inviting us to address separately different "syndrome patterns" solutions might have. For e.g., a quadruplet resulting as:

$$
\begin{gathered}
n_1 = X_0 X_1 X_4 X_7 \\
n_2 = X_0 X_1 \overline{X_4 X_7} \\
n_3 = X_0 X_1 X_4 \overline{X_7}\\
n_4 = X_0 \overline{X_1} X_4 \overline{X_7} \\
X_2=X_3=X_5=X_6=1 + i0
\end{gathered}
$$

Will be said to have a "${\{0,1,4,7\}}$" syndrome pattern. It is therefore natural to ask the following questions:

- What syndrome patterns exist in the known families of [3] or in sporadic solutions.
- Are there some syndrome patterns that we can rule out altogether?
- Conversely, are some patterns more promising then others?

One should also note that in the context of the $E_4$ problem, the integers are divided into a pair of pairs:

$$
\{\{n_1, n_2\},\{n_3, n_4\}\}
$$

And the $L_2$ condition is agnostic to the order within the pairs and between them. This imposes an equivalence on some of the syndrome patterns, inviting us to canonize the syndrome patterns with respect to this equivalence.

> **TODO (canonization lemma — makes the syndrome census well-defined):** Prove
> a unique canonical syndrome representative under unit twists, relabeling, and
> the pair/pair-of-pairs symmetry. Without it the "syndrome patterns" statistic
> in §7 is not well-defined. Also characterize and saturate away the
> **degenerate locus** of $W=0$ (repeated roots, coincident $n_i$, zero
> coordinates) before any zero-counting use. (`drafts/critical_review.md` §3.)

### The $L_{2,d}$ Condition and Fourth Powers

The $L_{2,d}$ term corresponds to the real part of the fourth powers of the Gaussian integers:

$$
L_{2,d}(a+ib, c+id) = \text{Re}\left\{(a + ib)^4 + (c + id)^4\right\}
$$

Applying the Brahmagupta-Fibonacci parameterization:

$$
\begin{gathered}
L_{2,d}(gh, g\overline{h}) = \text{Re}\left\{( gh )^4 + (g\overline{h})^4\right\}  = \text{Re}\left\{g^4(h^4 + \overline{h^4})\right\} = 2\text{Re}\left\{g^4\right\}\text{Re}\left\{h^4\right\} \\
= 2(m^2 - n^2 - 2mn)(m^2 - n^2 + 2mn)(p^2 - q^2 - 2pq)(p^2 - q^2 + 2pq)
\end{gathered}
$$

We believe that the factorization of the 8-th degree polynomial to four 2-nd degree components can be useful in both analytical and numerical studies of the search.

> **TODO ($\mathbb{Z}[\sqrt2]$ theorem — upgrades a §7 Finding from data to
> proof):** Each quartic factor $m^2 \pm 2mn - n^2 = (m\pm n)^2 - 2n^2$ is a
> value of the discriminant-$8$ norm form of $\mathbb{Z}[\sqrt2]$; hence every
> odd prime dividing $L_{2,d}$ is $\equiv \pm1 \pmod 8$. Add the smoothness
> explanation (forced factorization into four quadratic forms of size
> $\sim|L_{2,d}|^{1/4}$). (`drafts/critical_review.md` §4(iv).)

#### Single Constraint Formulation of $E_4$ - The $L_{2,d}$ Prism

The single-constraint formulation for $E_4$ becomes:

$$
\begin{gathered}
W(\bold{X}) = L_{2,d}(n_1(\bold{X}), n_2(\bold{X})) - L_{2,d}(n_3(\bold{X}), n_4(\bold{X})) = 0 \Rightarrow \\
\ \\
\text{Re}\left\{(X_0 X_1 X_2 X_3)^4\right\}\text{Re}\left\{(X_4 X_5 X_6 X_7)^4\right\} = \text{Re}\left\{(X_0 \overline{X_3} X_4 \overline{X_7})^4\right\}\text{Re}\left\{(X_1 \overline{X_2} X_5 \overline{X_6})^4\right\}
\end{gathered}
$$

This refactors the degree 32 equation into degree 8 terms, each involving only 8 variables, thus potentially simplifying numerical searches.

## 6. Algorithms

### Implementation-Aware Enumeration of $E_4$ Solutions

**Goal:** Find all $E_4$ solutions for a given norm $N$.

**Stage 1**

- Factorize $N$ and decompose $q_i \equiv 1 \pmod 4$ into Gaussian factors.
- The $p_i \equiv 3 \pmod 4$ primes and the factor of 2 only contribute global scaling or rotation so we omit them.

**Stage 2**

- Construct the set $U$ of all Gaussian integers with norm $\tilde{N} = q_1^{b_1} \dots q_s^{b_s}$. We compute $U$ in a tree-like fashion, branching over each factor and its conjugate.
  - Start with $U = \{1\}$
  - For each $q_1$ in $\tilde{N}$'s factorization (including multiplicities): $U \ \longleftarrow \ \mathrm{q_1}U  \bigcup \overline{\mathrm{q_1}}U$
- Canonize elements to the first octant ($a + ib: 0 \le b \le a$) and remove duplicates.
- Compute $r_{x,y} = (xy)^2$ (using two integer multiplications)

**Stage 3**

- Initialize a hash-table $H$.
- For each unordered pair $x+iy, u+iv \in U^2$
  - Compute key $k = r_{x,y} + r_{u,v}$
  - Try to insert $H[k]$, a collision yields a solution.

**Complexity and Implementation Details**

- Runtime and memory grow exponentially with the number of factors ($O(2^{2n})$).
- Large integer arithmetic can be optimized using reduced precision (e.g., $\mathbb{Z}/2^{32}\mathbb{Z}$) to detect collisions, followed by full-precision verification.
- In stage 2, We compute a total of $2^{n+1}$ gaussian integer multiplications, each such multiplication would normally cost 4 integer multiplications. We can utilize the fact that we also multiply with the conjugate to reduce this number to 2.
- $r_{x,y}$ was chosen with respect to $L_{2,b}$ and costs two multiplications. We could equivalently choose it to correspond to $L_{2,d}$ which would simply be the real component of the final integer, provided we precompute all factors to the fourth power (attractive in the reduced precision approach).
- The complexity of factorizing $N$ is not discussed since in our context it will rarely be an arbitrary large random number.
- Computing the decompositions of the gaussian factors can be done as a preliminary stage and reused.
- We are likely to enumerate large sets of norms, therefore we should not find it hard to parallelize the search. Furthermore, jointly enumerating certain norms together allows us to amortize some of the work, for e.g., if the norms share factors that saves some of arithmetic in Stage 2.

> **TODO (additional algorithms — see `drafts/critical_review.md` §6 and
> Front 1):**
>
> - **Non-squarefree extension** of the enumerator (the surjectivity gap in §5);
>   every exclusion claim depends on it.
> - **Top-collision** hash-join on $(L_1,L_2,L_3)$ across catalog and families.
> - **Bottom-lift** divisor join with the parity wheel ($\mu$ up to
>   $10^5\text{–}10^6$).
> - **Ladder congruence pre-sieve** (the forced divisors of §4) baked in
>   downstream.
> - **Symbolic family intersection** via resultants/Gröbner (blocked on
>   identifying [3]'s families symbolically).
> - **Hardening protocol** for the exclusion certificate: two independent
>   implementations, synthetic solution injection, the family self-test
>   (re-find [3]'s members incl. non-squarefree $k$), $\mathbb{F}_p$
>   density cross-checks, pre-registered bounds, checkpointed logs.

## 7. Computational Results — The Exhaustive Campaign

### Analysis of Known $E_4$ Solutions

We examined the following known solutions:

- The first family of [3] for $k=[1,30)$.
- The second family of [3] for $k=[5,60)$.
- Sporadic solutions 4 to 8 for table 1 of [3].
- A new solution found by our enumeration method.

For each solution we performed the following analysis:

- Canonized the solution:
  - Scaled to remove common factors.
  - Rotated to the first octant.
  - Transformed such that $L_1$ will be odd.
  - Sorted the solution pairs in descending order.
- Computed $L_1$ and $L_{2,d}$. If they were not too large, we factorized them.
- Computed the canonized syndrome pattern by factoring each of the four Gaussian integers and analyzed factorizations.

For e.g., our new solution exhibits the following properties:

$$
\begin{gathered}
A=252885 \ B=46703 \ C=195203 \ D=167415 \\
E=249703 \ F=61485 \ G=209985 \ H=148453 \\
L_1 = 33065996717 = 13^1\cdot 89^1\cdot 173^1\cdot 233^1\cdot 709^1 \\
L_{2,d}  = -912832201535971887688 = -1 \cdot 2^3\cdot 7^1\cdot 17^1\cdot 23^2\cdot 31^1\cdot 41^1\cdot 47^1\cdot 137^1\cdot 401^1\cdot 601^1\cdot 919^1 \\
\text{Syndrome Pattern} = \{0, 1, 3, 5\} \\
\end{gathered}
$$

We can draw the following insights:

- None of the 90 unique solutions we examined shared a value of $L_1$ or $L_{2,d}$, thus failing to pair up to create an $E_5$ solution.
- Unique Syndrome patterns (calculated for the smallest 18 solutions) are:

$$
\begin{gathered}
\{0, 1, 5, 6\}, \{0, 1, 3, 5, 6\}, \{0, 1, 4, 6\} \\
\{0, 1, 4, 5, 6\}, \{0, 1, 3, 4, 5, 6\}, \{0, 1, 2, 4, 5, 7\}, \\
\{0, 1, 2, 4\}, \{0, 1, 3, 4, 6\}, \{0, 1, 3, 5\}, \{0, 1, 3, 4, 5\} \\
\end{gathered}
$$

- The factorization of $L_1$ and $L_{2,d}$ appears to be surprisingly smooth (comparable to random integers of same size) and rarely contain multiplicities for odd primes.
  - Smoothness is probably related to the fact that under the parameterization, they are a product of 2 and 4 quadratic polynomials respectively.
- Furthermore, the primes composing $L_{2,d}$ are always $\equiv \pm 1 \pmod 8$ (In $L_1$ they are  $\equiv 1 \pmod 4$ which aligns with $L_1$ being a norm of a Gaussian integer).

> **TODO (the campaign — this section's headline experiment):** The 90-solution
> analysis above is a pilot, not a campaign. Fill in (see
> `drafts/critical_review.md` §7 and R1–R8):
>
> - Two **pre-registered regimes kept rigorously separate**: (i) exhaustive over
>   all canonical $L_1 \le B$ ($B$ from a pilot, order $10^{10}$–$10^{12}$),
>   which licenses exclusion; (ii) targeted champion norms (many small split
>   primes, $s \approx 20$ on GPU), which supports discovery but not
>   exhaustiveness.
> - The **explicit $E_5$ exclusion statement**: an exhaustive $E_4$ catalog to
>   $L_1 \le B$ excludes any $E_5$ with $L_1 \le B$, since both constituent
>   $E_4$'s share $L_1$.
> - Catalog + deduplication counts (mod the full symmetry group), the syndrome
>   census via the canonization algorithm, $(L_1, L_{2,d})$ collision statistics,
>   smoothness distributions against the forced-divisor baseline, density counts
>   against the §8 local model.
> - Cross-check per-norm counts against the $\mathbb{F}_p$ model to catch bugs.

## 8. A Density Heuristic

> **TODO (new section):** Multidegree / circle-method count (see
> `drafts/critical_review.md` §4(iii), §8). Start at $E_4$ (Tsai et al. [5] cover
> $E_3$):
>
> - $E_4$ is a type $(2,2,2,4)$ complete intersection in $\mathbb{P}^7$
>   (canonical class $\mathcal{O}(2)$, general type); circle-method exponent
>   $H^{8-10}=H^{-2}$ — yet infinite families exist, so the families live on
>   special rational subvarieties (Lang–Vojta).
> - $E_5$ is type $(2^7,4^3,8)$ in $\mathbb{P}^{15}$ ($\mathcal{O}(18)$, very
>   general type); exponent $H^{16-34}=H^{-18}$: essentially no accidental
>   solutions.
> - State as a conjecture with constants from the $\mathbb{F}_p$ run (§6/R7);
>   hedge on smoothness.

## 9. Future Directions

> **TODO:** Expand each into a paragraph; prioritize the ring pivot (highest
> expected value). Full drafts in `drafts/critical_review.md` "Future directions".

- **The ring pivot ($\mathbb{Z}[i]$).** The congruum step $L \pm u = \text{squares}$ parameterizes through factorizations since $x^2+y^2$ splits, so per-norm candidate sets grow like Gaussian divisor counts and level-2 collisions become far likelier; any $E_5$ over $\mathbb{Z}[i]$ would be the first ideal size-16 solution in any ring. Transplanting CMSV's rotation ansatz to size 16 over $\mathbb{Z}[i]$ is an independent second attack.
- **Geometry of the solution varieties.** Singular loci and low-degree rational curves on the $(2,2,2,4)$ complete intersection would mechanize family generation and ground §8; under Lang–Vojta, $E_5$ solutions off a special locus are finite.
- **Arithmetic dynamics.** $E$-towers are totally rational preimage trees for non-autonomous quadratic compositions $y^2 - L_i$; even an abc-conditional bound on maximal full-splitting depth would be a genuine theorem extending §3.
- **Completing the local theory.** A full 2-adic structure theorem and everywhere-local solubility at depth 5, upgrading §4 and §8.
- **Sharper required divisors.** Derive a constant $C^{E}_{16}$ strictly containing CMSV's $C'_{16}$ by exploiting correlations between sibling radicals — new arithmetic that would tighten §4 (blocked on reading CMSV [4] §3 for possible foreclosure).
- **Complexity packaging.** State precisely what a skew 32-gem does to the BCSS-variant landscape, with the geometry as a conditional explanation for why skew gems stop at $n=4$; extend the BMR [1] gem records.
- **Quaternion triage.** The Hurwitz-quaternion route to $L_{2,a}$ (a sum of four fourth-powers) faces the 24-element unit group and a non-multiplicatively-defined square-coordinate locus; one time-boxed week, then a remark either way.
- **CRT approach.** Solve modulo various primes and lift to integers.

## Appendix — Verified Exemplar

> **TODO:** Report the full recovered tower for the new solution above and
> verify it independently via the harness (§6 hardening). Per
> `drafts/critical_review.md` §2 the recovered values are
> $L_2 = 489628056848329146064$,
> $L_3 = 175480010455650701584492675662518592000000$,
> $\sqrt{L_4} = 40042900368028062136207226327668992000000$,
> with $\sqrt{L_4} = 2^{14}\cdot3^6\cdot5^6\cdot7^2\cdot11^2\cdot13\cdot17\cdot19
> \cdot29\cdot37^2\cdot43\cdot61\cdot107\cdot109\cdot139\cdot263\cdot389\cdot499$.
> Confirm the degree-16 polynomial identity $\prod(x^2-r_k^2) =
> (((x^2-L_1)^2-L_2)^2-L_3)^2-L_4$ and the ladder divisibilities before citing.

## References

> **TODO:** Replace bare URLs with full citations (authors, titles, years).
> Reference [3] must be identified from `sources/Bre08.pdf` (blocks the family
> intersection in §6). Add [4] CMSV and [5] Tsai et al. properly.

- [1] Borchert, McKenzie, Reinhardt, *Few Product Gates But Many Zeros* (MFCS 2009). `sources/BMR09.pdf`. <https://www.researchgate.net/publication/220975993_Few_Product_Gates_But_Many_Zeros>
- [2] Sum of Squares Function. <https://mathworld.wolfram.com/SumofSquaresFunction.html>
- [3] *(Experimental Mathematics, ~2009 — identify from `sources/Bre08.pdf`)*. <https://projecteuclid.org/journalArticle/Download?urlid=em%2F1243429952>
- [4] Coppersmith, Mossinghoff, Scheinerman, VanderKam (arXiv:2304.11254). `sources/CSMV23.pdf`. **TODO:** full citation.
- [5] Tsai et al. (2026), symmetric locus of the ideal degree-three PTE problem. `sources/TLT26.pdf`. **TODO:** full citation.
