# A Search for Normal 32-gem motivated by the PTE problem

- [A Search for Normal 32-gem motivated by the PTE problem](#a-search-for-normal-32-gem-motivated-by-the-pte-problem)
  - [Abstract](#abstract)
  - [Introduction](#introduction)
    - [The $E\_n$ Problem Statement](#the-e_n-problem-statement)
    - [Necessary and Sufficient Conditions for an $E\_n$ Solution](#necessary-and-sufficient-conditions-for-an-e_n-solution)
    - [Key Property - Solution Construction](#key-property---solution-construction)
  - [The Connection to Gaussian Integers](#the-connection-to-gaussian-integers)
    - [Factorization Analysis](#factorization-analysis)
  - [Joint Parameterization of $2^n$ Norm-like Gaussian Integers](#joint-parameterization-of-2n-norm-like-gaussian-integers)
    - [The $n=1$ Case and Application to the $E\_3$ Problem](#the-n1-case-and-application-to-the-e_3-problem)
    - [The $n=2$ Case and Application to the $E\_4$ Problem](#the-n2-case-and-application-to-the-e_4-problem)
    - [The General Case](#the-general-case)
  - [Formulation of the $E\_4$ Problem Using a Single Constraint](#formulation-of-the-e_4-problem-using-a-single-constraint)
  - [The $L\_{2,d}$ Condition and Fourth Powers](#the-l_2d-condition-and-fourth-powers)
    - [Single Constraint Formulation of $E\_4$ - The $L\_{2,d}$ Prism](#single-constraint-formulation-of-e_4---the-l_2d-prism)
  - [Implementation-Aware Enumeration of $E\_4$ Solutions](#implementation-aware-enumeration-of-e_4-solutions)
    - [Pinpointing a Particular $L\_1$](#pinpointing-a-particular-l_1)
    - [](#)
  - [More Ideas](#more-ideas)
  - [Appendix](#appendix)
  - [References](#references)

## Abstract

We investigate the search integers $(P,Q,R,S,T)$ such that the polynomial $((((x - P)^2 - Q)^2 - R)^2 - S)^2 - T^2$ factors into linear factors. Such a finding will provide an ideal symmetric solution for the Prouhet-Tarry-Escott (PTE) problem with degree $n=15$, extending the current state of the art currently at $n=11$. We propose several research directions which, to our knowledge are novel. We present both theoretical insights and an analysis of implementation details for extensive numerical enumeration.

- [ ] Our main strategy is:
- [ ] Our novel result is:

## Introduction

### The $E_n$ Problem Statement

Given integers $L_1, L_2, \dots, L_n$, we define the associated **$E_n$ Polynomial** recursively as:

$$E_1(x; L) := x^2 - L$$
$$E_n(x; L_1, \dots, L_n) := E_1(E_{n-1}(x; L_1, \dots, L_{n-1}); L_n)$$

The **$E_n$ Problem** is to find integers such that the polynomial $E_n(x; L_1, \dots, L_n)$ factors into linear factors over the integers. We refer to the tuple $(L_1, L_2, \dots, L_n)$ as an **$E_n$ solution**.

### Necessary and Sufficient Conditions for an $E_n$ Solution

The roots $r_1, \dots, r_{2^n}$ of an $E_n$ polynomial are arranged in a binary tree structure. Each node of the tree is associated with a constraint that must be satisfied. The first three constraints are presented below:

$$
\begin{gathered}
r_{2k}^2 + r_{2k+1}^2 = 2L_1 \\
(r_{4k}^2 - r_{4k+1}^2)^2 + (r_{4k+2}^2 - r_{4k+3}^2)^2 = 8L_2 \\
((r_{8k}^2 - r_{8k+1}^2)^2 - (r_{8k+2}^2 - r_{8k+3}^2)^2)^2 + ((r_{8k+4}^2 - r_{8k+5}^2)^2 - (r_{8k+6}^2 - r_{8k+7}^2)^2)^2 = 32L_3
\end{gathered}
$$

We refer to a constraint at the $n$-th level as the **$L_n$ condition**. Note that an $L_{k>1}$ condition possesses certain degrees of freedom arising from the lower conditions in the tree. In [1], an equivalent set of conditions is introduced under the name **Litter Conditions**. Useful augmentations of the $L_2$ condition include:

$$
\begin{gathered}
L_{2,a} := r_0^4 + r_1^4 + r_2^4 + r_3^4 = 4(L_1^2 + L_2)  \\
L_{2,b} := r_0^2r_1^2 + r_2^2r_3^2 = 2(L_1^2 - L_2)  \\
L_{2,c} := (r_0^2 - L_1)^2 + (r_2^2 - L_1)^2 = 2L_2  \\
L_{2,d} := (r_0^2 - r_1^2)^2 - 4r_0^2r_1^2 + (r_2^2 - r_3^2)^2 - 4r_2^2r_3^2 = 8(2L_2 - L_1^2)  \\
\end{gathered}
$$

### Key Property - Solution Construction

If two $E_n$ solutions share identical $L_1, \dots, L_{n-1}$ coefficients, they can be used to construct an $E_{n+1}$ solution. This suggests a strategy for finding $E_5$ solutions by first finding many $E_4$ solutions, and then find a pair of solutions taht share $P$, $Q$ and $R$.

## The Connection to Gaussian Integers

Since the $L_1$ condition for all root pairs implies:
$$
r_{2k}^2 + r_{2k+1}^2 = 2L_1 = 2L_1\\
$$

we can treat the root pairs as Gaussian integers of equal norm. When searching for roots that satisfy higher order $L_n$ conditions, we can confine our search to the complex circle of radius $2L_1$. By fixing $P$ and examining its factorization, we can:

1. Detenmine the number of points with norm $2L_1$ (namely the sum of two squares function, $r_2(n)$).
2. Efficiently ($O(n)$) enumerate all such points (see Appendix).
3. Parameterize sets of norm-like Gaussian integers even without a concrete $L_1$ at hand.

### Factorization Analysis

Let $n$ be a norm factoring as:

$$
n = 2^{a_0}p_1^{2a_1} \dots p_r^{2a_r}q_1^{b_1} \dots q_s^{b_s}
$$

- The $q_i$ are odd rational primes such that $q_i \equiv 1 \pmod 4$. They decompose as $q_i = \mathrm{q}_i\overline{\mathrm{q}}_i = (s_i + it_i)(s_i - it_i)$.
- The $p_i$ are odd rational primes such that $p_i \equiv 3 \pmod 4$. These remain inert in the Gaussian integers and must appear with even multiplicity. Their existence implies a common factor to all roots; thus, we generally assume the norm is free of such primes.
- The factor 2 contributes either a global scaling or a $45^\circ$ rotation, neither of which changes the essential structure of the roots.

Any Gaussian integer $x = a + ib$ with norm $a^2 + b^2 = n$ satisfies (up to units):

$$
x = \prod _{i=1}^s \mathrm{q}_i^{e_i} \overline{\mathrm{q}}_i^{b_i - e_i}
$$

This is characterized by a vector $\mathbf{e} = (e_1, \dots, e_s)$ where $0 \le e_i \le b_i$.

For the sake of simplicity, in the next sections we only treat norms without multiplicities. This allows us to view $\mathbf{e}$ as a binary vector. In that setting, taking the logical NOT of the $e_i$ is equivalent to taking the complex conjugate of $\mathrm{q}_i$.

## Joint Parameterization of $2^n$ Norm-like Gaussian Integers  

We show it is possible to parameterize any set of $2^n$ Gaussian integers sharing the same norm using $2^{2^n}$ variables. We are interested in essentially different solutions, therefore we pay special attention to the common factors and order of the roots.  

### The $n=1$ Case and Application to the $E_3$ Problem

Let $x, y$ be norm-like Gaussian integers associated with vectors $e_x$ and $e_y$. Let $g = m + in$ be the product of factors where $e_x, e_y$ are identical and, $h = p + iq$ where they differ.

$$
\begin{gathered}
x = gh = mp - nq + i(np - mq) \\
y = g\overline{h} = mp + nq + i(np + mq) \\
\end{gathered}
$$

This yields the Brahmagupta–Fibonacci parameterization of $a^2 + b^2 = c^2 + d^2$.Since the $L_1$ condition is the only constraint for the $E_3$ problem, every instance of this parameterization is a valid $E_3$ solution.

### The $n=2$ Case and Application to the $E_4$ Problem

To parameterize four norm-like Gaussian integers $n_1, n_2, n_3, n_4$, we identify $2^3 = 8$ distinct "syndromes" (subsets of indices where factors are conjugated relative to $n_1$). Let $X_i = s_i + it_i$ be the product of all factors matching the $i$-th syndrome. Following this observation, we can write down the four Gaussian integers as:

$$
\begin{gathered}
n_1 = X_0 X_1 X_2 X_3 X_4 X_5 X_6 X_7 \\
n_2 = X_0 X_1 X_2 X_3 \overline{X_4 X_5 X_6 X_7} \\
n_3 = X_0 X_1 \overline{X_2 X_3} X_4 X_5 \overline{X_6 X_7}\\
n_4 = X_0 \overline{X_1} X_2 \overline{X_3} X_4 \overline{X_5} X_6 \overline{X_7} \\
\end{gathered}
$$

Each of the coefficients $A,B,C,D,E,F,G,H$ can be thought of a degree 8 polynomial in 16 variables.

### The General Case

This generalizes to any number of integers (not only powers of two). The complexity of the parameterization makes it hard to utilize. For e.g., in light of the $E_5$ problem, we can parameterize a quadruplet of norm-like Gaussian integers. The expressions for the $L_2$ and $L_3$ conditions would result in polynomials with $2^8 = 256$ variables of degree no less than $512$ and $1024$ respectively.

## Formulation of the $E_4$ Problem Using a Single Constraint

The $L_2$ "Q condition" ($A^2B^2 + C^2D^2 = E^2F^2 + G^2H^2$) translates to finding the zeros of a degree 32 polynomial in 16 variables. This is the first known reduction of the $E_4$ problem to a single constraint. To the best of our knowledge, this is the first reduction of the $E_4$ problem into a single constraint problem (as opposed to the original formulation which involves 4 equations).

Investigating this polynomial is a difficult task due to it's size. The actual polynomial expression is far too big to be displayed here, and not even easily manipulated in a computer program. The following observation can be used to simplify the expression: Our quadruplet will not always contain all possible syndromes $X_0,  \dots  , X_7$, inviting us to address separately different "syndrome patterns" solutions might have. For e.g. a quadruplet resulting as:

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
  - [ ] Perform the analysis here or in the appendix.
- Are there some syndrome patterns that we can rule out altogether?
  - [ ] Discuss only a few examples.
- Conversely, are some patterns more promising then others?

One should also note that in the context of the $E_4$ problem, the integers are divided into a pair of pairs:

$$
\{\{n_1, n_2\},\{n_3, n_4\}\}
$$

And the Q condition is agnostic to order within the pairs and in between them. This imposes an equivalence on some of the syndrome patterns, inviting us to canonize the syndrome patterns w.r. to this equivalence.

## The $L_{2,d}$ Condition and Fourth Powers

The $L_{2,d}$ term corresponds to the real part of the fourth powers of the Gaussian integers:

$$
L_{2,d} = \text{Re}\left\{(a + ib)^4 + (c + id)^4\right\}
$$

Applying the $n=1$ parameterization:

$$
\begin{gathered}
L_{2,d}(gh, g\overline{h}) = \text{Re}\left\{( gh )^4 + (g\overline{h})^4\right\}  = \text{Re}\left\{g^4(h^4 + \overline{h^4})\right\} = 2\text{Re}\left\{g^4\right\}\text{Re}\left\{h^4\right\} \\
= 2(m^2 - n^2 - 2mn)(m^2 - n^2 + 2mn)(p^2 - q^2 - 2pq)(p^2 - q^2 + 2pq)
\end{gathered}
$$

We believe that the factorization of the 8-th degree polynomial to four 2-nd degree components can be useful in both in analytical and numerical study of the search.

### Single Constraint Formulation of $E_4$ - The $L_{2,d}$ Prism

The single-constraint formulation for $E_4$ becomes:

$$
\begin{gathered}
L_{2,d}(n_1, n_2) = L_{2,d}(n_3, n_4) \Rightarrow \\
\ \\
\text{Re}\left\{(X_0 X_1 X_2 X_3)^4\right\}\text{Re}\left\{(X_4 X_5 X_6 X_7)^4\right\} = \text{Re}\left\{(X_0 X_2 X_4 X_6)^4\right\}\text{Re}\left\{(X_1 X_3 X_5 X_7)^4\right\}
\end{gathered}
$$

This refactors the degree 32 equation into degree 8 terms, potentially simplifying numerical searches.

## Implementation-Aware Enumeration of $E_4$ Solutions

### Pinpointing a Particular $L_1$

**Goal:** Find all $E_4$ solutions for a given integer $P$.

**Stage 1**

- Factorize $P$ and decompose $q_i \equiv 1 \pmod 4$ into Gaussian factors.
- The $p_i \equiv 3 \pmod 4$ and the factor of 2 only contribute global scaling or rotation so we omit them in this context.

**Stage 2**

- Construct the set $U$ of all Gaussian integers with norm $\tilde{P} = q_1^{b_1} \dots q_s^{b_s}$. We compute $U$ in a tree-like fashion, branching over each factor and its conjugate.
  - Start with $U = \{1\}$
  - For each $q_1$ in $\tilde{P}$'s factorization (including multiplicity):
    - $U \ \longleftarrow \ \mathrm{q_1}U  \bigcup \overline{\mathrm{q_1}}U$
- Canonize elements to the first octant ($a + ib: 0 \le b \le a$) and remove duplicates.
- Compute $r_{x,y} = (xy)^2$ (using two integer multiplications)

**Stage 3**

- Initialize a hash-table $H$.
- For each unordered pair $x+iy, u+iv \in U^2$
  - Compute key $k = r_{x,y} + r_{u,v}$
  - Try to insert $H[k]$, a collision yields a solution.

**Complexity and Implementation Details**

- Runtime and memory grow exponentially with the number of factors ($O(2^{2n-1})$).
- Large integer arithmetic can be optimized using reduced precision (e.g., $\mathbb{Z}/2^{32}\mathbb{Z}$) to detect collisions, followed by full-precision verification.
- In stage 2, We compute a total of $O(2 \cdot 2^n)$ gaussian integer multiplications, each such multiplication, would normally cost 4 integer multiplications. We can utilize the fact that we also multiply with the conjugate to reduce this number to 2.
- $r_{x,y}$ was chosen with respect to $L_{2,b}$ and costs two multiplications. We could equivalently choose it to correspond to $L_{2,d}$ which would simply be the real component of the final integer, provided we precompute all factors to the fourth power (attractive in the reduced precision approach).
- Stage 3 complexity is quadratic
- The complexity of factorizing $P$ is not discussed since in our context it will rarely be an arbitrary large random number.
- Computing the decompositions of the gaussian factors can be done as a preliminary stage and reused.

###

## More Ideas

- Providing an estimate for the $E_{4/5}$ solution **density** to bound enumeration efforts.
- Systematically explore solution families via generalized linear constraints, extending the techniques of [3].
- Apply the Brahmagupta-Fibonacci parameterization to the $L_{2,b}$ or $L_{2,c}$ constraints, since they also can be viewed as sum of squares.
- The ring of Hurwitz quaternions admit a norm which is a sum of four squares. We can use it to parameterize the expression for $L_{2,a}$, which is a sum of four fourth-powers.
- Apply a Chinese Remainder Theorem (CRT) approach: solve modulo various primes and lift to integers.

## Appendix

- [ ] Delete if empty

## References

- [1] <https://www.researchgate.net/publication/220975993_Few_Product_Gates_But_Many_Zeros>
- [2] <https://mathworld.wolfram.com/SumofSquaresFunction.html>
- [3] <https://projecteuclid.org/journalArticle/Download?urlid=em%2F1243429952>
