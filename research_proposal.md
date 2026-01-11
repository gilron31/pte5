# A Search for Normal 32-gem motivated by the PTE problem

- [A Search for Normal 32-gem motivated by the PTE problem](#a-search-for-normal-32-gem-motivated-by-the-pte-problem)
  - [Abstract](#abstract)
  - [Introduction](#introduction)
    - [The $E\_n$ Problem Statement](#the-e_n-problem-statement)
    - [Necessary and Sufficient conditions for an $E\_n$ Solution](#necessary-and-sufficient-conditions-for-an-e_n-solution)
    - [Key Property - Solution Construction](#key-property---solution-construction)
  - [The Connection to Gaussian Integers](#the-connection-to-gaussian-integers)
  - [Joint Parameterization of $2^n$ norm-like Gaussian Integers](#joint-parameterization-of-2n-norm-like-gaussian-integers)
    - [The $n=1$ Case and Application to the $E\_3$ Problem](#the-n1-case-and-application-to-the-e_3-problem)
    - [The $n=2$ Case and Application to the $E\_4$ Problem](#the-n2-case-and-application-to-the-e_4-problem)
    - [The General Case](#the-general-case)
  - [Formulation of the $E\_4$ Problem using a single constraint](#formulation-of-the-e_4-problem-using-a-single-constraint)
    - [Ruling out Certain Syndrome patterns](#ruling-out-certain-syndrome-patterns)
    - [Characterization of Known $E\_4$ solution to syndromes](#characterization-of-known-e_4-solution-to-syndromes)
  - [The $L\_{2,d}$ Condition and Fourth Powers of Gaussian Integers](#the-l_2d-condition-and-fourth-powers-of-gaussian-integers)
    - [Single Constraint Formulation of $E\_4$ - The $L\_{2,d}$ Prism](#single-constraint-formulation-of-e_4---the-l_2d-prism)
  - [Implementation Aware Enumeration of $E\_4$ solutions](#implementation-aware-enumeration-of-e_4-solutions)
    - [Pin-Point a Particular $L\_1$](#pin-point-a-particular-l_1)
    - [](#)
  - [More Ideas](#more-ideas)
  - [Appendix](#appendix)
  - [References](#references)

## Abstract

We are interested in finding integers $(P,Q,R,S,T)$ such that the polynomial $((((x - P)^2 - Q)^2 - R)^2 - S)^2-T^2$ factors into linear factors. Such a finding will provide a solution for the PTE problem with $n=16$, extending the current state of the art currently at $n=12$. We propose several research directions which to the best of our knowledge are novel. We will go into both theoretical insights, as well as provide an analysis of possible implementation details whenever extensive enumeration is discussed.

- [ ] Our main strategy is:
- [ ] Our novel result is:

## Introduction

### The $E_n$ Problem Statement

Given integers $L_1, L_2,  \dots  , L_n$ we define the associated **$\bold{E_n}$ Polynomial** recursively as:

$$
E_1(x; L) := x^2 - L \\
E_n(x; L_1,  \dots , L_n) := E_1(E_{n-1}(x; L_1,  \dots , L_{n-1}); L_n) \\
$$

An **$\bold{E_n}$ Problem** is the problem of finding such integers so the polynomial $E_n(x; L_1,  \dots , L_n)$ factors into linear factors over the integers. We will therefore call $(L_1, L_2,  \dots  , L_n)$ an **$\bold{E_n}$ solution**.

### Necessary and Sufficient conditions for an $E_n$ Solution

One can show that the roots $r_1, \dots, r_{2^n}$  of an $E_n$ polynomial are arranged in a binary tree structure. Each node of the tree can be associated with a constraint that needs to be satisfied. An illustration of the first three constraints is written below:

$$
\begin{gathered}
r_{2k}^2 + r_{2k+1}^2 = 2L_1 \\
(r_{4k}^2 - r_{4k+1}^2)^2 + (r_{4k + 2}^2 - r_{4k+3}^2)^2 = 8L_2\\
((r_{8k}^2 - r_{8k+1}^2)^2 - (r_{8k + 2}^2 - r_{8k+3}^2)^2)^2 + ((r_{8k + 4}^2 - r_{8k+5}^2)^2 - (r_{8k + 6}^2 - r_{8k+7}^2)^2)^2 = 32L_3\\
\end{gathered}
$$

We will refer to a constraint an the n'th level as the $\bold{L_n}$ **condition**. One should note that an $L_{k>1}$ condition is enjoys certain degrees of freedom arising from the lower conditions in the tree. In [1], an equivalent set of conditions is introduced under the name **Litter Conditions**. Noteworthy augmentations of the $L_2$ condition include:

$$
\begin{gathered}
L_{2,a} := r_0^4 + r_1^4 + r_2^4 + r_3^4 = 4(L_1^2 + L_2)  \\
\ \\
L_{2,b} := r_0^2r_1^2 + r_2^2r_3^2 = 2(L_1^2 - L_2)  \\
\ \\
L_{2,c} := (r_0^2 - L_1)^2 + (r_2^2 - L_1)^2 = 2L_2  \\
\ \\
L_{2,d} := (r_0^2 - r_1^2)^2 - 4r_0^2r_1^2 + (r_2^2 - r_3^2)^2 - 4r_2^2r_3^2 = 8(2L_2 - L_1^2)  \\
\end{gathered}
$$

### Key Property - Solution Construction

If we find two $E_n$ solutions which have identical  $L_1,  \dots  , L_{n-1}$ coefficients, then we can construct a $E_{n+1}$ solution. We can apply this insight as a strategy for finding $E_5$ solutions by first finding many $E_4$ solutions, and then find a pair of solutions which share $P$, $Q$ and $R$.

## The Connection to Gaussian Integers

Since the $L_1$ condition for all root pairs implies:
$$
r_{2k}^2 + r_{2k+1}^2 = 2L_1 = 2P\\
$$

we can view the root pairs as a Gaussian integers, under the constraint that they are of equal norm. When searching for roots that satisfy higher order $L_n$ conditions, we can confine our search to the complex circle with a radius $2P$. By fixing $P$ and examining it's factorization, we can:

- [ ] Use conventional notation for this function.

1. Count the number of different points with norm $2P$ (namely the sum of two squares function).
2. Efficiently ($O(n)$) enumerate over all points (see Appendix).
3. Furthermore, even without a concrete $2P$ at hand, we can parameterize sets of norm-like Gaussian integers.

The following is a (rather straightforward) analysis of all gaussian integers of a particular norm.
Let's start by assuming that we have a set of gaussian numbers all with norm $n$ which factors as:

$$
n = 2^{a_0}p_1^{2a_1} \dots p_r^{2a_r}q_1^{b_1} \dots q_s^{b_s}
$$

- The $q_i$'s denote odd rational primes which decompose in the gaussian integers ($q_i = 1 \mod 4$). We shall denote their decomposition as $q_i = \mathrm{q}_i\overline{\mathrm{q}}_i = (s_i + it_i)(s_i - it_i)$, $0< t_i < s_i$.
- The $p_i$'s denote odd rational primes which remain primes in the gaussian integers ($p_i = 3 \mod 4$). They appear with even multiplicity because TODO. The existence of such factor in the norm directly implies a common factor to all roots, therefore in most cases we can assume w.l.o.g. that our norm is free of such primes.
- The factor $2$, contributes either a global scaling of all roots (if the power is even) or a simple mutation of the solutions (rotation by 45 degrees). In either cases it will not contribute an essential change to the roots so we will ignore it as well.

Therefore any gaussian integer $x = a + ib$ whose norm is $a^2 + b^2 = n$ must satisfy (up to units $\pm 1, \pm i$ which can be safely ignored):

$$
x = \prod _{i=1}^s \mathrm{q}_i^{e_i} \overline{\mathrm{q}}_i^{b_i - e_i}
$$

Therefore be characterized by a vector $\bold{e} = (e_1,  \dots ,  e_s)$ where $0 \le e_i \le b_i$.

For the sake of simplicity, in the next sections we choose to treat only the no-multiplicity case where $b_i = 1$ so we can view $\bold{e}$ as a binary vector. Is is worth noting that in that sense, taking the NOT of certain entries of that vector is equivalent of computing the conjugate of the corresponding $\mathrm{q}_i$. Treating the general case (with multiplicities) is equivalent to having repeated factors, when special care of ordering has to be taken to avoid degeneracies.

## Joint Parameterization of $2^n$ norm-like Gaussian Integers  

In this section we will show that it is possible to parameterize any set of $2^n$ Gaussian integers sharing the same norm using $2^{2^n}$ variables. Such parameterization can be helpful in finding $E_{n-2}$ solutions. It should be noted that we are interested in essentially different solutions, that means we pay special attention to the common factors and order of the roots.  

### The $n=1$ Case and Application to the $E_3$ Problem

Let $x = a + ib, y = c + id$ be two norm-like gaussian integers with norm $n$. Following the characterization described above, we can associate $x$ with a vector $e_x$ and $y$ with $e_y$. We can associate a gaussian integer $g = m + in$ with the entries where $e_x$ and $e_y$ are identical and a gaussian integer $h = p + iq$ where they differ. We can then write:

$$
\begin{gathered}
x = gh = mp - nq + i(np - mq) \\
y = g\overline{h} = mp + nq + i(np + mq) \\
\end{gathered}
$$

Which results in the known Brahmagupta–Fibonacci parameterization of the diophantine equation $a^2 + b^2 = c^2 + d^2$. [3] Completes the discussion in the context of $E_3$ problems. It turns out however, that the $L_1$ condition is the only constraints required for the fulfillment of the $E_3$ problem, therefore every instance of the parameterization is a valid $E_3$ solution.  

### The $n=2$ Case and Application to the $E_4$ Problem

We now consider the case where where we have four norm-like gaussian integers $n_1, n_2, n_3, n_4$ whose parameterization we seek. Same as before, we associate each with a binary vector $e_1, e_2, e_3, e_4$. The task of finding a concise parameterization now translates into finding disjoint subsets of indices, s.t. for each subset has a distinct "syndrome" over each of the four vectors. Each subset has a corresponds to a "base form" which is a gaussian integer, defined by the product of gaussian decomposable primes composing $e_1$ at the subset's indices. The "syndrome" of a given subset refers to for which of $e_2, e_3, e_4$ do the base appears conjugated. There are exactly 8 different syndromes, which we will denote by the numbers 0-7. Therefore we require exactly 8 subsets to fully parameterize our four integers. We will denote the resulting gaussian integers corresponding to each subset (and syndrome) with $X_i = s_i + i t_i$. We can therefore construct the following parameterization for 4 norm-like gaussian integers:

$$
\begin{gathered}
n_1 = A + iB = X_0 X_1 X_2 X_3 X_4 X_5 X_6 X_7 \\
n_2 = C + iD = X_0 X_1 X_2 X_3 \overline{X_4 X_5 X_6 X_7} \\
n_3 = E + iF = X_0 X_1 \overline{X_2 X_3} X_4 X_5 \overline{X_6 X_7}\\
n_4 = G + iH = X_0 \overline{X_1} X_2 \overline{X_3} X_4 \overline{X_5} X_6 \overline{X_7} \\
\end{gathered}
$$

Each of the coefficients $A,B,C,D,E,F,G,H$ of the resulting $n_i$ is a degree 8 polynomial in 16 variables.

### The General Case

While we won't go into the details here, it is easy to see how we can generalize our parameterization to any number of norm-like gaussian integers (not only powers of two). Specifically we could (in principle) write down the $L_2$ and $L_3$ conditions for the eight gaussian integers required for the satisfaction of the $E_5$ problem (it would require a staggering $2^8 = 256$ variables and the resulting constraints would be of degree no less than $2^7 \cdot 4 = 512$).

## Formulation of the $E_4$ Problem using a single constraint

Given the parameterization of the previous section one can be easily convinced that the task of satisfying the $L_2$ condition, the "Q condition":

$$
A^2B^2 + C^2D^2 = E^2F^2 + G^2H^2
$$

Now translates to finding the zeros of a degree 32 polynomial in 16 variables. To the best of our knowledge, this is the first reduction of the $E_4$ problem into a single constraint problem (as opposed to the original formulation which involves 4 equations).

Investigating this polynomial is a difficult task due to it's size. The actual polynomial expression is far too big to be displayed here, and not even easily manipulated in a computer program. We need ways to simplify this result which would still be useful, hence the following observation: our quadruplet will not always contain all possible syndromes $X_0,  \dots  , X_7$, so we can discuss the different "syndrome patterns" different solutions might have. For e.g. a quadruplet resulting as:

$$
\begin{gathered}
n_1 = X_0 X_1 X_4 X_7 \\
n_2 = X_0 X_1 \overline{X_4 X_7} \\
n_3 = X_0 X_1 X_4 \overline{X_7}\\
n_4 = X_0 \overline{X_1} X_4 \overline{X_7} \\
X_2=X_3=X_5=X_6=1 + i0
\end{gathered}
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

## The $L_{2,d}$ Condition and Fourth Powers of Gaussian Integers

Viewing the root pairs as Gaussian integers, the $L_{2,d}$ term can be viewed as the sum of the real part of the fourth powers of the two Gaussian integers:

$$
L_{2,d} = Re\left\{(a + ib)^4 + (c + id)^4\right\}
$$

We can take a step further and apply the Brahmagupta–Fibonacci parameterization here:

$$
\begin{gathered}
L_{2,d}(g= m + in,h = p + iq) = Re\left\{( gh )^4 + (g\overline{h})^4\right\}  = Re\left\{g^4(h^4 + \overline{h^4})\right\} = 2Re\left\{g^4\right\}Re\left\{h^4\right\} \\
= 2(m^2 - n^2 - 2mn)(m^2 - n^2 + 2mn)(p^2 - q^2 - 2pq)(p^2 - q^2 + 2pq)
\end{gathered}
$$

We believe that the factorization of the 8-th degree polynomial to four 2-nd degree components might be useful in both in analytical and numerical study of the search.

### Single Constraint Formulation of $E_4$ - The $L_{2,d}$ Prism

We can now apply equations TODO to rewrite the single-constraint parameterization of the $E_4$ problem as:

$$
\begin{gathered}
n_1 = A + iB = X_0 X_1 X_2 X_3 X_4 X_5 X_6 X_7 \\
n_2 = C + iD = X_0 X_1 X_2 X_3 \overline{X_4 X_5 X_6 X_7} \\
n_3 = E + iF = X_0 X_1 \overline{X_2 X_3} X_4 X_5 \overline{X_6 X_7}\\
n_4 = G + iH = X_0 \overline{X_1} X_2 \overline{X_3} X_4 \overline{X_5} X_6 \overline{X_7} \end{gathered}
$$

$$
\begin{gathered}
L_{2,d}(n_1, n_2) = L_{2,d}(n_3, n_4) \Rightarrow \\
\ \\
Re\left\{(X_0 X_1 X_2 X_3)^4\right\}Re\left\{(X_4 X_5 X_6 X_7)^4\right\} = Re\left\{(X_0 X_2 X_4 X_6)^4\right\}Re\left\{(X_1 X_3 X_5 X_7)^4\right\}
\end{gathered}
$$

Thus we have refactored our original equation of degree 32 in 16 variables, to smaller degree 8 terms in 8 variables. While encouraging, it is not immediately clear how this might be useful.

## Implementation Aware Enumeration of $E_4$ solutions

### Pin-Point a Particular $L_1$

- Inputs - An integer $P$
- Outputs - All $E_4$ solutions with $L_1 = P$
- Stage 1:
  - Factorize $P$ over the integers. Using the notation of section TODO $P = 2^{a_0}p_1^{a_1} \dots p_r^{a_r}q_1^{b_1} \dots q_s^{b_s}$. W.l.o.g. we can ignore the $2$ and the $p_i$ and re-introduce them in a post-process stage.
  - Decompose all $q_i$ to their gaussian factors $\mathrm{q_i}$.
- Stage 2:
  - Construct the set $U$ all Gaussian integers with norm $\tilde{P} = q_1^{b_1} \dots q_s^{b_s}$.
    - Start with $U = \{1\}$
    - For each $q_1$ in $\tilde{P}$'s factorization (including multiplicity):
      - $U \ \longleftarrow \ \mathrm{q_1}U  \bigcup \overline{\mathrm{q_1}}U$
  - It is important to canonize elements to the first eighth of the complex plane ($a + ib: 0 \le b \le a$) and then to remove duplicates.
  - For each $x + iy \in U$:
    - Compute $r_{x,y} = (xy)^2$ (using two integer multiplications)
- Stage 3:
  - Initialize a hash-table $H$.
  - For each unordered pair $x+iy, u+iv \in U^2$
    - Compute key $k = r_{x,y} + r_{u,v}$
    - Try to insert $H[k]$
      - If no previous result is stored, insert value $(x, y, u ,v)$.
      - If a previous result is stored, yield a collision $((x,y,u,v), H[k])$.

- Complexity analysis and Implementation details:
  - We omit the discussion about the factorization of $P$ since in our context it will rarely be an arbitrary large random number.
  - Computing the decompositions of the gaussian factors can be done as a preliminary stage and reused.
  - Since collisions are rare, but the integer arithmetic is heavy (especially multiprecision arithmetic), we can perform stages 2-3 with reduced precision, for e.g. $\mathbb{Z} /  2^{32}\mathbb{Z}$ or some other hardware-friendly ring. We might encounter false collisions that way but those would be easier to rule out in a secondary verification stage.
  - The runtime (and memory) complexity of this algorithm grows exponentially in $n$, the number of factors of $\tilde{P}$.
  - In stage 2, We compute a total of $O(2 \cdot 2^n)$ gaussian integer multiplications, each such multiplication, would normally cost 4 integer multiplications but we can reduce that number by a factor of 2 since we jointly compute the multiplication of an arbitrary complex number $a$ with another complex number $b$ and his conjugate $\overline{b}$. Computing $r_{x,y}$ costs another two multiplications per gaussian integer.
  - In stage 3, the compute cost and memory demands for the hash-table are both proportional to the number of unordered pairs $\sim 2^{2n-1}$.

###

## More Ideas

- Providing an estimate for the $E_{4/5}$ solution **density** can be a powerful tool in estimating the extent of enumeration needed to find them.
- More systematic exploration of solution families following/extending Bremner's techniques. For e.g.:
  - Imposing generalized linear constraints on the roots instead of the specialized (2-3) and (2-4) constraints.
- Parameterizing $L_2$ condition as a gaussian integer.
  - Both in standard form
  - and in "centralized form"
- Using Hurwitz quaternions for the $L_2$ condition.
- Gaining possible insight the extension tower $\mathbb{Q}(\zeta_4) \supset \mathbb{Q}(i) \supset \mathbb{Q}$ In the sense that $\mathbb{Q}(\zeta_8)$ can be viewed as manipulating pair of gaussian integers.
- Applying a "CRT" approach to the problem as a whole. Solving modulu different primes and then construct an integer solution using CRT.
- Trying to disprove the existance of $E_5$ solutions.

## Appendix

- [ ] Delete if empty

## References

- [1] <https://www.researchgate.net/publication/220975993_Few_Product_Gates_But_Many_Zeros>
- [2] <https://mathworld.wolfram.com/SumofSquaresFunction.html>
- [3] <https://projecteuclid.org/journalArticle/Download?urlid=em%2F1243429952>
