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

For integers $L_1, \dots, L_n$ let $E_n(x) = ((\cdots(x^2 - L_1)^2 - L_2)^2 \cdots)^2 - L_n$, a polynomial of degree $2^n$; we call $(L_1, \dots, L_n)$ an **$E_n$ solution** when $E_n$ splits completely over $\mathbb{Z}$. Such polynomials are exactly the *skew $2^n$-gems* of Borchert, McKenzie and Reinhardt, and — in pairs sharing $L_1, \dots, L_{n-1}$ — they yield ideal symmetric solutions of the Prouhet–Tarry–Escott (PTE) problem of size $2^{n-1}$. No $E_5$ solution (PTE size $16$) is known, and none exists below the extensive symmetric searches of Coppersmith, Mossinghoff, Scheinerman and VanderKam. We first correct a folklore correspondence: $E_n$ solutions are precisely the **sym-perfect** pairs of [\[BMR\]](#ref-BMR), a *proper* subclass of symmetric ideal solutions — Tarry's 1913 degree-$7$ solution is symmetric ideal but is not an $E_4$. As the headline result we prove the first explicit arithmetic lower bounds for an $E_5$: a **divisibility ladder** cascades the PTE-constant divisors of [\[CMSV\]](#ref-CMSV) through the nested radicals of the tower, forcing any $E_5$ over $\mathbb{Z}$ to satisfy $\sqrt{L_5} \ge C'_{16}/2 \approx 2.44 \times 10^{27}$, hence $L_5 \ge 5.9 \times 10^{54}$; combined with [\[CMSV\]](#ref-CMSV)'s exhaustive size-$16$ search, its largest root exceeds $850$ and $L_1 > 3.6 \times 10^{5}$. Supporting this we give a self-similar recursion identifying $E_n$ solutions with congruum-lifted $E_{n-1}$ value systems; a multiplicative parameterization over $\mathbb{Z}[i]$ reducing the $E_4$ conditions to a single polynomial constraint $W(\mathbf X) = 0$ (making explicit the scalar constraint of Dilcher and the variety of Bremner); and a proof that the fourth-power invariant $L_{2,d}$ factors through the norm form of $\mathbb{Z}[\sqrt 2]$, so its odd prime divisors are $\equiv \pm 1 \pmod 8$. We then give per-norm enumeration and top-/bottom-collision search algorithms and report an exhaustive catalog of $E_4$ solutions.

## 1. Introduction

### The Prouhet–Tarry–Escott problem

The **Prouhet–Tarry–Escott (PTE) problem** asks for two distinct multisets of integers $A = \{a_1, \dots, a_m\}$ and $B = \{b_1, \dots, b_m\}$ of equal size with equal power sums up to some degree $k$:

$$
\sum_{i=1}^m a_i^\ell = \sum_{i=1}^m b_i^\ell, \qquad \ell = 1, 2, \dots, k.
$$

Equivalently, $\prod_i (x - a_i) - \prod_i (x - b_i)$ has degree at most $m - k - 1$. Since a nontrivial solution forces $m \ge k+1$, the extremal case $k = m - 1$ is called **ideal**; a solution is **symmetric** if each set is stable under negation about a common center. The problem is classical — it traces to Euler and Goldbach and was studied by Prouhet (1851), Tarry and Escott (1910s); see [\[RN\]](#ref-RN) and [\[BoIn\]](#ref-BoIn) for surveys.

Ideal solutions are difficult to produce and the state of the art is narrow. Over $\mathbb{Z}$, ideal solutions are known for every size $m \le 10$ and for $m = 12$, but **none is known for $m = 11$ or for any $m \ge 13$** [\[CMSV\]](#ref-CMSV), [\[RN\]](#ref-RN); the largest known is size $12$ (degree $11$). Recently Coppersmith, Mossinghoff, Scheinerman and VanderKam [\[CMSV\]](#ref-CMSV) developed the divisibility theory of the *PTE constant* and conducted exhaustive searches for symmetric ideal solutions: they found nothing at sizes $13$–$16$ (the size-$16$ search reaching height $850$), confirmed exactly two primitive symmetric solutions at size $12$, and discovered new ideal solutions over $\mathbb{Z}[i]$ and other imaginary-quadratic rings at sizes up to $12$ — but nothing above size $12$ in any ring.

### The iterated-squaring ansatz

This paper studies a rigid subclass of symmetric ideal solutions arising from *iterated squaring*. For integers $L_1, \dots, L_n$, define the **$E_n$ polynomial** by $E_1(x) = x^2 - L_1$ and $E_k(x) = E_{k-1}(x)^2 - L_k$; explicitly

$$
E_n(x) = \big( \cdots ((x^2 - L_1)^2 - L_2)^2 \cdots \big)^2 - L_n,
$$

a monic polynomial of degree $2^n$. We call $(L_1, \dots, L_n)$ an **$E_n$ solution** when $E_n$ factors completely into linear factors over $\mathbb{Z}$. The $E_3$ case was parameterized by Dilcher [\[Dil\]](#ref-Dil), who reduced it to writing an integer as a sum of two squares in two ways; the same reduction underlies the recent density results of Tsai, Lee and Takahashi [\[TLT\]](#ref-TLT) for degree-$3$ symmetric PTE (our $E_3$). The $E_4$ case is exactly Bremner's problem of splitting $(((X^2 - P)^2 - Q)^2 - R)^2 - S^2$ [\[Bre\]](#ref-Bre); he constructed two infinite families via elliptic curves and posed the existence of $E_5$ (five nested squares) as **open**. Notably, Dilcher had earlier *claimed $E_4$ impossible* — an assertion refuted by Bremner and by an example of Symes.

Iterated-squaring polynomials also carry algebraic-complexity content. In the terminology of Borchert, McKenzie and Reinhardt [\[BMR\]](#ref-BMR), an $E_n$ solution is a **skew $2^n$-gem**: a $\{+,-,\times\}$-circuit computing a degree-$2^n$ polynomial with $2^n$ distinct integer roots using the minimum $n$ product gates. Skew $2^n$-gems for infinitely many $n$ would refute the Bürgisser $L$-conjecture and the Blum–Cucker–Shub–Smale $\tau$-conjecture (hence $\mathrm{P}_{\mathbb{C}} \ne \mathrm{NP}_{\mathbb{C}}$), and connect to integer factoring. Gems are known through $n = 4$; a skew $32$-gem ($E_5$) would be the first at $n = 5$, and [\[BMR\]](#ref-BMR) report failing to find one.

### Contributions

We keep the focus on the $E_n$ ansatz itself, rather than on the size-$16$ record. After fixing definitions and the PTE correspondence (§2), our contributions are, **bounds first**:

1. **The first explicit lower bounds for an $E_5$ solution (§4).** A divisibility ladder cascades [\[CMSV\]](#ref-CMSV)'s required divisors of the PTE constant through the nested radicals of the tower, forcing $\sqrt{L_5} \ge C'_{16}/2 \approx 2.44 \times 10^{27}$ and $L_5 \ge 5.9 \times 10^{54}$ for any $E_5$ over $\mathbb{Z}$; together with [\[CMSV\]](#ref-CMSV)'s exhaustive size-$16$ search this gives largest root $> 850$ and $L_1 > 3.6 \times 10^5$.
2. **Corrected structure theory (§2–§3).** $E_n$ solutions are exactly the *sym-perfect* subclass of [\[BMR\]](#ref-BMR), properly contained in symmetric ideal solutions; and a self-similar recursion (§3) identifies an $E_n$ solution with an $E_{n-1}$ solution whose entries are *congrua* of a common $L_1$, yielding a second, bottom-up search axis.
3. **A single-constraint parameterization (§5).** A multiplicative parameterization of equal-norm Gaussian integers turns the $E_4$ conditions into a single polynomial identity $W(\mathbf X) = 0$, making explicit Dilcher's scalar constraint and Bremner's variety; and the fourth-power invariant $L_{2,d}$ is shown to factor through the norm form of $\mathbb{Z}[\sqrt 2]$, forcing its odd prime divisors to be $\equiv \pm 1 \pmod 8$.
4. **Algorithms and an exhaustive $E_4$ catalog (§6–§7)**, with a density heuristic (§8) and directions (§9).

All numeric claims in §2–§5 are reproduced by `drafts/verify_claims.py`. Throughout we take care to attribute the ansatz, the tree conditions, the correspondence, and the single-constraint reduction to their sources; the genuinely new ingredients are the ladder bounds, the multiplicative parameterization, the $\mathbb{Z}[\sqrt2]$ theorem, and the enumeration campaign.

## 2. The $E_n$ Problem and its PTE Correspondence

### The root tree and the tree conditions

Writing $E_n(x; L_1, \dots, L_n)$ when we wish to display the parameters, recall $E_1(x; L) = x^2 - L$ and $E_n(x; L_1, \dots, L_n) = E_1\!\big(E_{n-1}(x; L_1, \dots, L_{n-1}); L_n\big)$. Our primary goal is a **nontrivial $E_5$ solution**. We first record the root structure of a split $E_n$ and then its equivalence with symmetric ideal PTE solutions.

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
((r_{8k}^2 - r_{8k+1}^2)^2 - (r_{8k+2}^2 - r_{8k+3}^2)^2)^2 \\
\qquad {} + ((r_{8k+4}^2 - r_{8k+5}^2)^2 - (r_{8k+6}^2 - r_{8k+7}^2)^2)^2 = 32L_3
\end{gathered}
$$

We refer to a constraint at the $n$-th level as the **$L_n$ condition**. Note that an $L_{n>1}$ condition possesses certain degrees of freedom arising from the lower conditions in the tree. These are exactly the **litter conditions** of [\[BMR, Def. 5.2\]](#ref-BMR) (whence the name): on the squared roots $r_0^2, \dots, r_{2^{n-1}-1}^2$, labelling each internal node of the tree by the product of the leaves below it, the litter conditions require the two children of every node at a fixed level to have equal sum; the $L_1$ and $L_2$ rows above are [\[BMR\]](#ref-BMR)'s first two litter conditions verbatim. Useful augmentations of the $L_2$ condition include:

$$
\begin{gathered}
L_{2,a} := r_0^4 + r_1^4 + r_2^4 + r_3^4 = 4(L_1^2 + L_2)  \\
L_{2,b} := r_0^2r_1^2 + r_2^2r_3^2 = 2(L_1^2 - L_2)  \\
L_{2,c} := (r_0^2 - L_1)^2 + (r_2^2 - L_1)^2 = 2L_2  \\
L_{2,d} := (r_0^2 - r_1^2)^2 - 4r_0^2r_1^2 + (r_2^2 - r_3^2)^2 - 4r_2^2r_3^2 = 8(2L_2 - L_1^2)  \\
\end{gathered}
$$

The four forms $L_{2,a},\dots,L_{2,d}$ are equivalent to the $L_2$ condition modulo the $L_1$ condition; we single out $L_{2,d}$ in §5.

### The PTE correspondence

The equivalence between $E_n$ solutions and PTE solutions is subtler than sometimes stated, and we make it precise. Suppose $E_n$ splits completely, and consider its top-level factorization $E_n = (E_{n-1} - \sqrt{L_n})(E_{n-1} + \sqrt{L_n})$. Let $S$ and $T$ be the root multisets of the two degree-$2^{n-1}$ factors $E_{n-1}(x; L_1, \dots, L_{n-1} \mp \sqrt{L_n})$. Because $\big(\prod_{s \in S}(x-s)\big) - \big(\prod_{t\in T}(x-t)\big) = -2\sqrt{L_n}$ has degree $0$, the pair $(S,T)$ is an **ideal** PTE solution of size $2^{n-1}$; and since $E_{n-1}$ is even, $S$ and $T$ are symmetric about $0$. This is [\[BMR, Cor. 3.5\]](#ref-BMR): **every $E_n$ solution yields a symmetric ideal PTE solution of size $2^{n-1}$.** Equivalently, an $E_n$ solution *is* a pair of $E_{n-1}$ solutions sharing $L_1, \dots, L_{n-2}$ (the two factors above), which is the point of view of the Key Property below.

The converse fails: not every symmetric ideal solution comes from a tower. [\[BMR, Thm. 3.7\]](#ref-BMR) characterizes those that do as the **sym-perfect** pairs — a recursively defined subclass — and shows the inclusion is proper. Concretely, the obstruction is already the $L_1$ condition. Consider the symmetric ideal degree-$7$ solution obtained by centering Tarry's 1913 solution $\{0,4,9,23,27,41,46,50\}$, $\{1,2,11,20,30,39,48,49\}$: its positive halves are $\{2,16,21,25\}$ and $\{5,14,23,24\}$. For this to be an $E_4$, each half would have to split into root pairs of equal norm $2L_1$; but the pairwise norm sums are $2^2+25^2 = 629$ and $16^2+21^2 = 697$ — unequal — so the $L_1$ condition fails and **this solution is not an $E_4$** (already observed in [\[BMR, p. 7\]](#ref-BMR)). Hence $E_n$ solutions form a *proper* subclass of symmetric ideal solutions; the tower structure is a genuine restriction, and it is exactly what makes the problem parameterizable (§5) while implying that null results for $E_5$ do not settle the symmetric size-$16$ problem, nor conversely.

### Key Property — solution construction

Dually to the top split, **if two $E_n$ solutions share $L_1, \dots, L_{n-1}$, their union is an $E_{n+1}$ solution**: writing them as $E_n(x; L_1, \dots, L_{n-1}, M_1)$ and $E_n(x; L_1, \dots, L_{n-1}, M_2)$, one checks that $E_{n-1}(x; L_1, \dots, L_{n-1})$ is the common inner tower and $E_{n+1} = E_n^2 - L_{n+1}$ splits with $\sqrt{L_{n+1}} = (M_2 - M_1)/2$ and new coefficient $L_n = (M_1 + M_2)/2$. This gives the natural top-down strategy for $E_5$: find many $E_4$ solutions and seek a pair sharing $L_1, L_2, L_3$. §3 gives a complementary bottom-up strategy.

## 3. The Tower Recursion

The tree conditions of §2 have a self-similar structure: the conditions on an $E_n$ root system, read one level up, are the conditions of an $E_{n-1}$ root system. This is the recursive content of [\[BMR\]](#ref-BMR)'s sym-perfect condition, which we make explicit in a form convenient for the algorithms of §6 and the bounds of §4. Order the $2^{n-1}$ positive roots $r_0, \dots, r_{2^{n-1}-1}$ tree-wise, so that $(r_{2k}, r_{2k+1})$ is a bottom-level pair with $r_{2k}^2 + r_{2k+1}^2 = 2L_1$, and define the **derived sequence**

$$
u_k := \frac{r_{2k}^2 - r_{2k+1}^2}{2}, \qquad k = 0, \dots, 2^{n-2}-1 .
$$

**Proposition 3.1 (self-similarity).** *Let $(r_k)$ be the tree-ordered positive roots of an $E_n$ solution with parameters $(L_1, \dots, L_n)$, and let $(u_k)$ be its derived sequence. Then:*

*(a) (congruum) for each $k$, both $L_1 + u_k = r_{2k}^2$ and $L_1 - u_k = r_{2k+1}^2$ are perfect squares;*

*(b) (self-similarity) the values $(|u_k|)$ are the tree-ordered positive roots of an $E_{n-1}$ solution with parameters $(L_2, L_3, \dots, L_n)$.*

*Proof sketch.* Part (a) is immediate from $r_{2k}^2 = L_1 \pm u_k$. For (b), substitute $r_{2k}^2 - r_{2k+1}^2 = 2u_k$ into the $L_2$ condition of §2: $(r_{4k}^2 - r_{4k+1}^2)^2 + (r_{4k+2}^2 - r_{4k+3}^2)^2 = 8L_2$ becomes $(2u_{2k})^2 + (2u_{2k+1})^2 = 8L_2$, i.e. $u_{2k}^2 + u_{2k+1}^2 = 2L_2$ — precisely the $L_1$ condition for $(u_k)$ with parameter $L_2$. The same substitution carries the $L_{j+1}$ condition on $(r_k)$ to the $L_j$ condition on $(u_k)$ for every $j$, shifting parameters $(L_2, \dots, L_n)$ into place. $\square$

We verified Proposition 3.1 for $n = 4$ on the exemplar of §7. Its derived sequence

$$
u = (30884826508,\ 5038214492,\ 29285591492,\ 11027703508)
$$

satisfies $u_0^2 + u_1^2 = u_2^2 + u_3^2 = 2L_2$, and its squared values solve $\prod_k (x^2 - u_k^2) = ((x^2 - L_2)^2 - L_3)^2 - L_4$ — an $E_3$ solution with parameters $(L_2, L_3, L_4)$ — while $L_1 \pm u_k$ are perfect squares for all $k$ (`drafts/verify_claims.py`). A careful accounting of the powers of $2$ in the higher $L_j$-normalizations, and the general induction, are routine and deferred to the full write-up.

> **TODO:** Give the general $2$-power bookkeeping and the parity corollary
> ($u_k \equiv 0 \pmod 4$ under an odd-$L_1$ normalization); a clean statement
> needs the normalization lemma of §5.

**Two search axes.** Part (a) says each $u_k$ is a *congruum* of $L_1$ — a value with $L_1 \pm u_k$ both square. This yields a bottom-up dual to the Key Property:

**Corollary 3.2 (bottom-lift criterion).** *An $E_{n+1}$ solution lies above a given $E_n$ solution with positive roots $(r'_k)$ if and only if there is a scaling $\mu$ and an integer $L$ such that $L + \mu r'_k$ and $L - \mu r'_k$ are both perfect squares for every $k$.* Indeed, writing $L + \mu r'_k = a_k^2$ and $L - \mu r'_k = b_k^2$ gives $a_k^2 - b_k^2 = 2\mu r'_k$, so each condition is a divisor enumeration of $2\mu r'_k$ subject to the shared value $L = (a_k^2 + b_k^2)/2$; screening a whole $E_4$ catalog for liftability is therefore an $8$-way hash/divisor join. Whereas the Key Property splits a tower at the *top* (pairing two $E_n$ that agree on $L_1, \dots, L_{n-1}$), Corollary 3.2 builds from the *bottom*; together with Proposition 3.1 they expose every intermediate cut point, enabling a meet-in-the-middle search (§6).

## 4. The Divisibility Ladder and $E_5$ Lower Bounds

This is the headline of the paper: an $E_5$ solution, if it exists, is forced to be arithmetically enormous. The mechanism combines the tower structure of §2 with the divisibility theory of the PTE constant developed by [\[CMSV\]](#ref-CMSV).

### The PTE constant and its required divisors

For an ideal PTE solution $(A, B)$ of size $n$, the polynomial difference $\prod_i(x - a_i) - \prod_i(x - b_i)$ is a constant, denoted $C_n(A, B)$. Let $C_n$ (resp. $C'_n$) be the greatest common divisor of $C_n(A, B)$ over all ideal (resp. symmetric ideal) solutions of size $n$. Kleiman, Rees–Smyth, Caley, Filaseta–Markovich and [\[CMSV\]](#ref-CMSV) established that these constants are highly composite; [\[CMSV, Table 1\]](#ref-CMSV) records the required divisors. The two rows we need are

$$
C'_8 \supseteq 2^6 \cdot 3^3 \cdot 5^2 \cdot 7^2 \cdot 11 \cdot 13, \qquad
C'_{16} \supseteq 2^{11} \cdot 3^6 \cdot 5^4 \cdot 7^3 \cdot 11^2 \cdot 13^2 \cdot 17 \cdot 19 \cdot 23 \cdot 29 \cdot 37 \cdot 41 \cdot 43 \cdot 53,
$$

the second combining the $C_{16}$ divisors with the additional symmetric-case primes $29, 37, 41, 43, 53$. Also relevant is $C'_4 = 2^2 \cdot 3^2 = 36$.

### The ladder

**Theorem 4.1 (divisibility ladder).** *Let $(L_1, \dots, L_n)$ be an $E_n$ solution. At each internal node of the tower, the two child subtowers form a symmetric ideal PTE solution whose constant equals twice the corresponding nested radical; consequently $C'_{2^{k}}/2$ divides the level-$k$ radical for every node. In particular, for the top split, $C'_{2^{n-1}} \mid 2\sqrt{L_n}$.*

*Proof.* Fix a node whose subtower is an $E_{k+1}$ with parameters $(M_1, \dots, M_{k+1})$; its top split (§2) is a symmetric ideal size-$2^{k}$ solution $(S, T)$ with $\prod_{S}(x - s) - \prod_{T}(x - t) = -2\sqrt{M_{k+1}}$. Thus $C_{2^{k}}(S, T) = -2\sqrt{M_{k+1}}$, and since $(S,T)$ is symmetric, $C'_{2^{k}} \mid 2\sqrt{M_{k+1}}$; as $C'_{2^k}$ is even this gives $C'_{2^{k}}/2 \mid \sqrt{M_{k+1}}$, and $\sqrt{M_{k+1}}$ is the radical at that node. $\square$

We confirmed Theorem 4.1 on the $E_4$ exemplar of §7 (`drafts/verify_claims.py`): $C'_8/2 = 151{,}351{,}200 \mid \sqrt{L_4}$, and $C'_4/2 = 18 \mid \sqrt{L_3 \pm \sqrt{L_4}}$ for both signs (both radicals being perfect squares, as the tree requires).

### Lower bounds for $E_5$

**Corollary 4.2 (first explicit $E_5$ bounds).** *Any $E_5$ solution over $\mathbb{Z}$ satisfies*

$$
\begin{gathered}
\sqrt{L_5} \;\ge\; \frac{C'_{16}}{2} = 2\,437\,428\,918\,743\,498\,865\,144\,960\,000 \approx 2.44 \times 10^{27}, \\
L_5 \ge (C'_{16}/2)^2 \approx 5.9 \times 10^{54}.
\end{gathered}
$$

*Proof.* Theorem 4.1 at the top node gives $C'_{16} \mid 2\sqrt{L_5}$, so $2\sqrt{L_5}$ is a positive multiple of $C'_{16}$, whence $\sqrt{L_5} \ge C'_{16}/2$. $\square$

This bound on the tower depth is complemented by a bound on the tower *height* coming from [\[CMSV\]](#ref-CMSV)'s exhaustive search. Their symmetric size-$16$ search found no solution of height $\le 850$. Since the top split of an $E_5$ is a symmetric ideal size-$16$ solution whose height is the largest root $r_{\max}$, we get $r_{\max} > 850$; and as $r_{\max}^2 = L_1 + \sqrt{L_2 + \cdots} \in [L_1, 2L_1]$ (the smallest root being real), this gives $L_1 > 850^2/2 = 361{,}250$. We record this separately because the two arguments are independent: the ladder bounds $L_5$ from *below* via local arithmetic; the exhaustive search bounds the height. (The ladder alone yields only $L_1 \ge (C'_{16}/2)^{1/8} \approx 2.7 \times 10^3$ through the chain $\sqrt{L_5} \le L_1^{8}$, weaker than the search.)

**A correction.** An earlier estimate of this bound used $C'_{16}/2 \approx 1.43 \times 10^{26}$; that value omits the prime $17$ (equivalently, it computes $C'_{16}/(2 \cdot 17)$). Since $17 \mid C_{16}$ in [\[CMSV, Table 1\]](#ref-CMSV), the correct value is the one above, larger by a factor of $17$. The exponents in $C'_{16}$ should be re-checked directly against [\[CMSV, Table 1\]](#ref-CMSV) before publication; the primes are what dominate the magnitude.

**Congruence sieve.** Beyond the magnitude bound, Theorem 4.1 forces specific prime content at every level of any $E_5$ tower (e.g. $18 \mid \sqrt{L_3 \pm \sqrt{L_4}}$ at each level-$3$ node). These forced divisors are a congruence pre-sieve for the searches of §6.

## 5. Gaussian Parameterization and the Single-Constraint Reduction

Since the $L_1$ condition for all root pairs implies $r_{2k}^2 + r_{2k+1}^2 = 2L_1$, we can treat the root pairs as Gaussian integers of equal norm. When searching for roots that satisfy higher order $L_n$ conditions, we can confine our search to the complex circle of radius $2L_1$. By fixing $L_1$ and examining its factorization, we can:

1. Determine the number of points with norm $2L_1$ (namely the sum-of-two-squares function $r_2(2L_1)$; its mean square drives the density heuristic of §8 [\[BC\]](#ref-BC)).
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

This yields the Brahmagupta–Fibonacci parameterization of $a^2 + b^2 = c^2 + d^2$. Since the $L_1$ condition is the only constraint for the $E_3$ problem, every instance of this parameterization is a valid $E_3$ solution, and the converse also holds — every $E_3$ solution is parameterized by some $g, h$. This recovers Dilcher's characterization of $E_3$ as an integer expressible as a sum of two squares in two ways [\[Dil\]](#ref-Dil), and the sum-of-two-squares reduction underlying the density count of [\[TLT\]](#ref-TLT) for degree-$3$ symmetric PTE.

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

Each of the coefficients $A,B,C,D,E,F,G,H$ can be thought of as a degree 8 polynomial in 16 variables $\mathbf{X} = \{s_i + it_i\}_{i=0}^{7}$.

#### The General Case

This generalizes to any number of integers. The exponential ($2^{n-1}$) growth of the number of possible syndromes rapidly makes the parameterizations hard to utilize. For e.g., in light of the $E_5$ problem, we can parameterize an octet of norm-like Gaussian integers. The expressions for the $L_2$ and $L_3$ conditions would result in dense polynomials with $2^8 = 256$ variables of degree no less than $512$ and $1024$ respectively.

### Formulation of the $E_4$ Problem Using a Single Constraint

The $L_2$ condition can be written as (Here is the $L_{2,b}$ form, but the result is agnostic to this augmentation):

$$
\begin{gathered}
W(\mathbf{X}) := A(\mathbf{X})^2 B(\mathbf{X})^2 +  C(\mathbf{X})^2 D(\mathbf{X})^2 - E(\mathbf{X})^2 F(\mathbf{X})^2 - G(\mathbf{X})^2 H(\mathbf{X})^2 = 0\\
\end{gathered}
$$

By the construction of the parameterization, the $L_1$ condition:

$$
\begin{gathered}
A(\mathbf{X})^2 + B(\mathbf{X})^2 = C(\mathbf{X})^2 + D(\mathbf{X})^2 \\
= E(\mathbf{X})^2 + F(\mathbf{X})^2 = G(\mathbf{X})^2 + H(\mathbf{X})^2\\
\end{gathered}
$$

is automatically satisfied. Therefore the $E_4$ problem now translates to finding the zeros of $W(\mathbf{X})$, a degree 32 polynomial in 16 variables.

That $E_4$ reduces to a *single* scalar constraint is not new: it is Dilcher's condition $a_1^2 b_1^2 + a_2^2 b_2^2 = a_3^2 b_3^2 + a_4^2 b_4^2$ among four equal-norm representations [\[Dil, (3-2)\]](#ref-Dil) (the same relation he used, incorrectly, to argue $E_4$ impossible), and the underlying locus is Bremner's degree-$32$ threefold in $\mathbb{P}^7$ cut by three quadrics and one quartic [\[Bre, §2\]](#ref-Bre). Our contribution is the *multiplicative* parameterization: by construction it trivializes the three $L_1$ (quadric) conditions, turning the single constraint into a polynomial identity $W(\mathbf X) = 0$ in the free variables $\mathbf X$, rather than a relation among constrained representations. This is what makes the syndrome analysis below and the enumeration of §6 possible.

Investigating $W(\mathbf{X})$ is a difficult task due to its size. The actual polynomial expression is far too big to be displayed here, and not even easily manipulated with a computer program. The following observation can be used to simplify the expression: Our quadruplet will not always contain all possible syndromes $X_0,  \dots  , X_7$, inviting us to address separately different "syndrome patterns" solutions might have. For e.g., a quadruplet resulting as:

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

- What syndrome patterns exist in the known families of [\[Bre\]](#ref-Bre) or in sporadic solutions.
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

The factorization of this degree-$8$ form into four quadratic factors has arithmetic consequences, which explain an empirical regularity of §7.

**Theorem 5.1 ($\mathbb{Z}[\sqrt2]$ constraint on $L_{2,d}$).** *For a primitive $E_4$ solution, every odd prime dividing $L_{2,d}$ is $\equiv \pm 1 \pmod 8$.*

*Proof.* Each quadratic factor above is $m^2 \pm 2mn - n^2 = (m \pm n)^2 - 2n^2$, a value $N_{\mathbb{Z}[\sqrt2]}(x + y\sqrt2) = x^2 - 2y^2$ of the norm form of $\mathbb{Z}[\sqrt2]$ (discriminant $8$), with $x = m \pm n$, $y = n$; likewise for the $p, q$ factors. Let $\pi$ be an odd prime dividing $L_{2,d} = 2\,\mathrm{Re}(g^4)\,\mathrm{Re}(h^4)$, hence dividing one such value $x^2 - 2y^2$. If $\pi \nmid y$, then $2 \equiv (x y^{-1})^2 \pmod \pi$ is a quadratic residue, so $\pi \equiv \pm 1 \pmod 8$. If $\pi \mid y$ then $\pi \mid x$, forcing a common factor of the corresponding Gaussian integer $g$ (or $h$); for a primitive solution this is excluded. $\square$

This upgrades to a theorem the empirical observation (§7) that the odd primes of $L_{2,d}$ are always $\equiv \pm1 \pmod 8$. It also explains the observed *smoothness*: $L_{2,d}$ is, up to the factor $2$, a product of four norm-form values each of size $\sim |L_{2,d}|^{1/4}$, so its prime factorization behaves like that of four independent integers a quarter the bit-length, which are correspondingly smoother than a single random integer of the same size. (The same mechanism, one level down, is Dilcher's observation that the $L_2$-difference is divisible by $2^7 \cdot 3^2$ [\[Dil\]](#ref-Dil).)

#### Single Constraint Formulation of $E_4$ - The $L_{2,d}$ Prism

The single-constraint formulation for $E_4$ becomes:

$$
\begin{gathered}
W(\mathbf{X}) = L_{2,d}(n_1(\mathbf{X}), n_2(\mathbf{X})) - L_{2,d}(n_3(\mathbf{X}), n_4(\mathbf{X})) = 0 \Rightarrow \\
\ \\
\text{Re}\left\{(X_0 X_1 X_2 X_3)^4\right\}\text{Re}\left\{(X_4 X_5 X_6 X_7)^4\right\} \\
\qquad = \text{Re}\left\{(X_0 \overline{X_3} X_4 \overline{X_7})^4\right\}\text{Re}\left\{(X_1 \overline{X_2} X_5 \overline{X_6})^4\right\}
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
>   identifying [\[Bre\]](#ref-Bre)'s families symbolically — now available, see `drafts/sources_insight.md` §1).
> - **Hardening protocol** for the exclusion certificate: two independent
>   implementations, synthetic solution injection, the family self-test
>   (re-find [\[Bre\]](#ref-Bre)'s members incl. non-squarefree $k$), $\mathbb{F}_p$
>   density cross-checks, pre-registered bounds, checkpointed logs.

## 7. Computational Results — The Exhaustive Campaign

### Analysis of Known $E_4$ Solutions

We examined the following known solutions:

- The first family of [\[Bre\]](#ref-Bre) for $k=[1,30)$.
- The second family of [\[Bre\]](#ref-Bre) for $k=[5,60)$.
- Sporadic solutions 4 to 8 for table 1 of [\[Bre\]](#ref-Bre).
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
L_{2,d}  = -912832201535971887688 \\
\qquad = -\,2^3\cdot 7\cdot 17\cdot 23^2\cdot 31\cdot 41\cdot 47\cdot 137\cdot 401\cdot 601\cdot 919 \\
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
> `drafts/critical_review.md` §4(iii), §8). Start at $E_4$ ([\[TLT\]](#ref-TLT) cover
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
- **Sharper required divisors.** Derive a constant $C^{E}_{16}$ strictly containing [\[CMSV\]](#ref-CMSV)'s $C'_{16}$ by exploiting correlations between sibling radicals — new arithmetic that would tighten §4. Their §3 constructs local solutions obstructing some divisibilities (any prime $p>n$ with $p\equiv\pm1\bmod n$ cannot be forced), so a genuine improvement must exploit the *towered* structure specifically; the 2-adic machinery of [\[FM\]](#ref-FM) is the natural starting point.
- **Complexity packaging.** State precisely what a skew 32-gem does to the BCSS-variant landscape, with the geometry as a conditional explanation for why skew gems stop at $n=4$; extend the [\[BMR\]](#ref-BMR) gem records.
- **Quaternion triage.** The Hurwitz-quaternion route to $L_{2,a}$ (a sum of four fourth-powers) faces the 24-element unit group and a non-multiplicatively-defined square-coordinate locus; one time-boxed week, then a remark either way.
- **CRT approach.** Solve modulo various primes and lift to integers.

## Appendix — Verified Exemplar

> **TODO:** Report the full recovered tower for the new solution above and
> verify it independently via the harness (§6 hardening). Per
> `drafts/critical_review.md` §2 the recovered values are
> $L_2 = 489628056848329146064$,
> $L_3 = 175480010455650701584492675662518592000000$,
> $\sqrt{L_4} = 40042900368028062136207226327668992000000$,
> whose factorization $\sqrt{L_4} = 2^{14}\cdot3^6\cdot5^6\cdot7^2\cdot11^2\cdot13\cdots499$
> (smooth, largest prime $499$) is listed in full in `drafts/verify_claims.py`.
> Confirm the degree-16 polynomial identity $\prod(x^2-r_k^2) =
> (((x^2-L_1)^2-L_2)^2-L_3)^2-L_4$ and the ladder divisibilities before citing.

## References

- []{#ref-BMR}**[BMR]** B. Borchert, P. McKenzie, K. Reinhardt. *Few product gates but many zeros.* MFCS 2009, LNCS 5734, Springer, 480–491; *Chicago J. Theoret. Comput. Sci.* 2013, art. 2, 1–22. (`sources/BMR09.pdf`)
- []{#ref-Bre}**[Bre]** A. Bremner. *When can $(((X^2-P)^2-Q)^2-R)^2-S^2$ split into linear factors?* *Experimental Mathematics* 17:4 (2008), 385–390. (`sources/Bre08.pdf`)
- []{#ref-BC}**[BC]** J. M. Borwein, K.-K. S. Choi. *On Dirichlet series for sums of squares.* *Ramanujan J.* 7 (2003), 95–127. (`sources/borwein-sums1.pdf`)
- []{#ref-BoIn}**[BoIn]** P. Borwein, C. Ingalls. *The Prouhet–Tarry–Escott problem revisited.* *L'Enseignement Math.* 40 (1994), 3–27.
- []{#ref-CMSV}**[CMSV]** D. Coppersmith, M. J. Mossinghoff, D. Scheinerman, J. M. VanderKam. *Ideal solutions in the Prouhet–Tarry–Escott problem.* arXiv:2304.11254 (2023). (`sources/CSMV23.pdf`)
- []{#ref-Dil}**[Dil]** K. Dilcher. *Nested squares and evaluations of integer products.* *Experimental Mathematics* 9:3 (2000), 369–372. (`sources/Dil00.pdf`)
- []{#ref-FM}**[FM]** M. Filaseta, M. Markovich. *Newton polygons and the Prouhet–Tarry–Escott problem.* *J. Number Theory* 174 (2017), 384–400. (`sources/PTEProblemPaper2016Submitted.pdf`)
- []{#ref-RN}**[RN]** S. Raghavendran, V. Narayanan. *The Prouhet–Tarry–Escott problem: A review.* *Mathematics* 7(3):227 (2019). (`sources/RN19.pdf`)
- []{#ref-TLT}**[TLT]** Y.-D. Tsai, J. Lee, F. Takahashi. *Arithmetic symmetry in ideal Prouhet–Tarry–Escott solutions.* arXiv:2606.07735 (2026). (`sources/TLT26.pdf`)
- []{#ref-Wro}**[Wro]** J. Wróblewski. *A collection of numerical solutions of multigrade equations related to the Prouhet–Tarry–Escott problem*, v7 (2009). (`sources/PTE007.pdf`)

*Still to add (cited within, PDFs not yet in `sources/`): Dorwart–Brown (1937), Chernick (1937), Rees–Smyth (1990), Kleiman (1975), Caley (2012/13). See `drafts/sources_insight.md` §7 for full details.*
