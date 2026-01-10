# A Search for Normal 32-gem motivated by the PTE problem

- [A Search for Normal 32-gem motivated by the PTE problem](#a-search-for-normal-32-gem-motivated-by-the-pte-problem)
  - [Abstract](#abstract)
  - [Introduction](#introduction)
    - [The $E\_n$ Problem Statement](#the-e_n-problem-statement)
    - [Necessary and Sufficient conditions for an $E\_n$ Solution - The Litter Conditions](#necessary-and-sufficient-conditions-for-an-e_n-solution---the-litter-conditions)
    - [Key Property - Solution Construction](#key-property---solution-construction)
  - [The Connection to Gaussian Integers](#the-connection-to-gaussian-integers)
  - [Joint Parameterization of $2^n$ norm-like Gaussian Integers](#joint-parameterization-of-2n-norm-like-gaussian-integers)
    - [The $n=1$ Case and Application to the $E\_3$ Problem](#the-n1-case-and-application-to-the-e_3-problem)
    - [The $n=2$ Case and Application to the $E\_4$ Problem](#the-n2-case-and-application-to-the-e_4-problem)
    - [The General Case](#the-general-case)
  - [Formulation of the $E\_4$ Problem using a single constraint](#formulation-of-the-e_4-problem-using-a-single-constraint)
    - [Ruling out Certain Syndrome patterns](#ruling-out-certain-syndrome-patterns)
    - [Characterization of Known $E\_4$ solution to syndromes](#characterization-of-known-e_4-solution-to-syndromes)
  - [Implementation Aware Enumeration of $E\_4$ solutions](#implementation-aware-enumeration-of-e_4-solutions)
  - [More Ideas](#more-ideas)
  - [Appendix](#appendix)
  - [References](#references)

## Abstract

We are interested in finding integers $(P,Q,R,S,T)$ such that the polynomial $((((x - P)^2 - Q)^2 - R)^2 - S)^2-T^2$ factors into linear factors. Such a finding will provide a solution for the PTE problem with $n=16$, extending the current state of the art currently at $n=12$. We propose several research directions which to the best of our knowledge are novel. We will go into both theoretical insights, as well as provide an analysis of possible implementation details whenever extensive enumeration is discussed.

- [ ] Our main strategy is:
- [ ] Our novel result is:

## Introduction

### The $E_n$ Problem Statement

Given integers $L_1, L_2, ... , L_n$ we define the associated **$\bold{E_n}$ Polynomial** recursively as:

$$
E_1(x; L) := x^2 - L \\
E_n(x; L_1, ..., L_n) := E_1(E_{n-1}(x; L_1, ..., L_{n-1}); L_n) \\
$$

An **$\bold{E_n}$ Problem** is the problem of finding such integers so the polynomial $E_n(x; L_1, ..., L_n)$ factors into linear factors over the integers. We will therefore call $(L_1, L_2, ... , L_n)$ an **$\bold{E_n}$ solution**.

### Necessary and Sufficient conditions for an $E_n$ Solution - The Litter Conditions

- [ ] Define roots.
- [ ] Define **K-th litter condition**
- [ ] Cover Bremner.
- [ ] Cover the other paper.

### Key Property - Solution Construction

If we find two $E_n$ solutions which have identical  $L_1, ... , L_{n-1}$ coefficients, then we can construct a $E_{n+1}$ solution. We can apply this insight as a strategy for finding $E_5$ solutions by first finding many $E_4$ solutions, and then find a pair of solutions which share $P$, $Q$ and $R$.

## The Connection to Gaussian Integers

Since the 1-st Litter condition for all root pairs implies:
$$
r_{2k}^2 + r_{2k+1}^2 = 2L_1 = 2P\\
$$

we can view the root pairs as a Gaussian integers, under the constraint that they are of equal norm. When searching for roots that satisfy higher order Litter conditions, we can confine our search to the complex circle with a radius $2P$. By fixing $P$ and examining it's factorization, we can:

- [ ] Use conventional notation for this function.

1. Count the number of different points with norm $2P$ (namely the sum of two squares function).
2. Efficiently ($O(n)$) enumerate over all points (see Appendix).
3. Furthermore, even without a concrete $2P$ at hand, we can parameterize sets of norm-like Gaussian integers.

The following is a (rather straightforward) analysis of all gaussian integers of a particular norm.
Let's start by assuming that we have a set of gaussian numbers all with norm $n$ which factors as:

$$
n = 2^{a_0}p_1^{2a_1}...p_r^{2a_r}q_1^{b_1}...q_s^{b_s}
$$

- The $q_i$'s denote odd rational primes which decompose in the gaussian integers ($q_i = 1 \mod 4$). We shall denote their decomposition as $q_i = \mathrm{q}_i\overline{\mathrm{q}}_i = (s_i + it_i)(s_i - it_i)$, $0< t_i < s_i$.
- The $p_i$'s denote odd rational primes which remain primes in the gaussian integers ($p_i = 3 \mod 4$). They appear with even multiplicity because TODO. The existence of such factor in the norm directly implies a common factor to all roots, therefore in most cases we can assume w.l.o.g. that our norm is free of such primes.
- The factor $2$, contributes either a global scaling of all roots (if the power is even) or a simple mutation of the solutions (rotation by 45 degrees). In either cases it will not contribute an essential change to the roots so we will ignore it as well.

Therefore any gaussian integer $x = a + ib$ whose norm is $a^2 + b^2 = n$ must satisfy (up to units $\pm 1, \pm i$ which can be safely ignored):

$$
x = \prod _{i=1}^s \mathrm{q}_i^{e_i} \overline{\mathrm{q}}_i^{b_i - e_i}
$$

Therefore be characterized by a vector $\bold{e} = (e_1, ...,  e_s)$ where $0 \le e_i \le b_i$.

For the sake of simplicity, in the next sections we choose to treat only the no-multiplicity case where $b_i = 1$ so we can view $\bold{e}$ as a binary vector. Is is worth noting that in that sense, taking the NOT of certain entries of that vector is equivalent of computing the conjugate of the corresponding $\mathrm{q}_i$. Treating the general case (with multiplicities) is equivalent to having repeated factors, when special care of ordering has to be taken to avoid degeneracies.

## Joint Parameterization of $2^n$ norm-like Gaussian Integers  

In this section we will show that it is possible to parameterize any set of $2^n$ Gaussian integers sharing the same norm using $2^{2^n}$ variables. Such parameterization can be helpful in finding $E_{n-2}$ solutions. It should be noted that we are interested in essentially different solutions, that means we pay special attention to the common factors and order of the roots.  

### The $n=1$ Case and Application to the $E_3$ Problem

Let $x = a + ib, y = c + id$ be two norm-like gaussian integers with norm $n$. Following the characterization described above, we can associate $x$ with a vector $e_x$ and $y$ with $e_y$. We can associate a gaussian integer $g = m + in$ with the entries where $e_x$ and $e_y$ are identical and a gaussian integer $h = p + iq$ where they differ. We can then write:

$$
x = gh = mp - nq + i(np - mq) \\
y = g\overline{h} = mp + nq + i(np + mq) \\
$$

Which results in the known parameterization of the diophantine equation $a^2 + b^2 = c^2 + d^2$. [3] Completes the discussion in the context of $E_3$ problems. It turns out however, that the 1-st Litter condition is the only constraints required for the fulfillment of the $E_3$ problem, therefore every instance of the parameterization is a valid $E_3$ solution.  

### The $n=2$ Case and Application to the $E_4$ Problem

We now consider the case where where we have four norm-like gaussian integers $n_1, n_2, n_3, n_4$ whose parameterization we seek. Same as before, we associate each with a binary vector $e_1, e_2, e_3, e_4$. The task of finding a concise parameterization now translates into finding disjoint subsets of indices, s.t. for each subset has a distinct "syndrome" over each of the four vectors. Each subset has a corresponds to a "base form" which is a gaussian integer, defined by the product of gaussian decomposable primes composing $e_1$ at the subset's indices. The "syndrome" of a given subset refers to for which of $e_2, e_3, e_4$ do the base appears conjugated. There are exactly 8 different syndromes, which we will denote by the numbers 0-7. Therefore we require exactly 8 subsets to fully parameterize our four integers. We will denote the resulting gaussian integers corresponding to each subset (and syndrome) with $X_i = s_i + i t_i$. We can therefore construct the following parameterization for 4 norm-like gaussian integers:

$$
n_1 = A + iB = X_0 X_1 X_2 X_3 X_4 X_5 X_6 X_7 \\
n_2 = C + iD = X_0 X_1 X_2 X_3 \overline{X_4 X_5 X_6 X_7} \\
n_3 = E + iF = X_0 X_1 \overline{X_2 X_3} X_4 X_5 \overline{X_6 X_7}\\
n_4 = G + iH = X_0 \overline{X_1} X_2 \overline{X_3} X_4 \overline{X_5} X_6 \overline{X_7} \\
$$

Each of the coefficients $A,B,C,D,E,F,G,H$ of the resulting $n_i$ is a degree 8 polynomial in 16 variables.

### The General Case

While we won't go into the details here, it is easy to see how we can generalize our parameterization to any number of norm-like gaussian integers (not only powers of two). Specifically we could (in principle) write down the 2-nd and 3-rd Litter conditions for the eight gaussian integers required for the satisfaction of the $E_5$ problem (it would require a staggering $2^8 = 256$ variables and the resulting constraints would be of degree no less than $2^7 \cdot 4 = 512$).

## Formulation of the $E_4$ Problem using a single constraint

Given the parameterization of the previous section one can be easily convinced that the task of satisfying the 2-nd Litter condition, the "Q condition":

$$
A^2B^2 + C^2D^2 = E^2F^2 + G^2H^2
$$

Now translates to finding the zeros of a degree 32 polynomial in 16 variables. To the best of our knowledge, this is the first reduction of the $E_4$ problem into a single constraint problem (as opposed to the original formulation which involves 4 equations).

Investigating this polynomial is a difficult task due to it's size. The actual polynomial expression is far too big to be displayed here, and not even easily manipulated in a computer program. We need ways to simplify this result which would still be useful, hence the following observation: our quadruplet will not always contain all possible syndromes $X_0, ... , X_7$, so we can discuss the different "syndrome patterns" different solutions might have. For e.g. a quadruplet resulting as:

$$
n_1 = X_0 X_1 X_4 X_7 \\
n_2 = X_0 X_1 \overline{X_4 X_7} \\
n_3 = X_0 X_1 X_4 \overline{X_7}\\
n_4 = X_0 \overline{X_1} X_4 \overline{X_7} \\
\\
\ \\
X_2=X_3=X_5=X_6=1 + i0
$$

Will be said to have a "${\{0,1,4,7\}}$" syndrome pattern. It is therefore natural to ask two questions:

- What syndrome patterns can we identify in the $E_4$ solutions we have discovered so far (either via enumeration or via the families of [3])?
- Are there some syndrome patterns that we can rule out altogether? are some more promising search grounds then others?

One should also note that in the context of the $E_4$ problem, the integers are divided into a pair of pairs:

$$
\{\{n_1, n_2\},\{n_3, n_4\}\}
$$

And the Q condition is agnostic to order within the pairs and in between them. This imposes an equivalence on some of the syndrome patterns.  

### Ruling out Certain Syndrome patterns

- [ ] Do

### Characterization of Known $E_4$ solution to syndromes

- [ ] Bremner families
- [ ] Sporadic solutions

## Implementation Aware Enumeration of $E_4$ solutions

- [ ] Write about:
  - The tree like $O(n)$ evaluation
  - meet in the middle
  - mult save tricks
  - precision lowering, two stage elimination

## More Ideas

- Providing an estimate for the $E_{4/5}$ solution **density** can be a powerful tool in estimating the extent of enumeration needed to find them.
- Imposing generalized linear constraints on the roots and searching for solution families with techniques similar to [3].
- Parameterizing the 2-nd Litter condition as a gaussian integer.
  - Both in standard form
  - and in "centralized form"
- Using Hurwitz quaternions for the 2-nd Litter condition.
- Gaining possible insight the extension tower $\mathbb{Q}(\zeta_4) \supset \mathbb{Q}(i) \supset \mathbb{Q}$ In the sense that $\mathbb{Q}(\zeta_4)$ can be viewed as manipulating pair of gaussian integers.

## Appendix

- [ ] Delete if empty

## References

- [1] <https://www.researchgate.net/publication/220975993_Few_Product_Gates_But_Many_Zeros>
- [2] <https://mathworld.wolfram.com/SumofSquaresFunction.html>
- [3] <https://projecteuclid.org/journalArticle/Download?urlid=em%2F1243429952>
