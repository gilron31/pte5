# A Search for a Normal 32-gem motivated by the PTE problem

## Abstract

We are interested in finding integers $(P,Q,R,S,T)$ such that the polynomial $((((x - P)^2 - Q)^2 - R)^2 - S)^2-T^2$ factors into linear factors. Such a finding will provide an example of the PTE conjecture with $n=15$, extending the current state of the art currently at $n=11$. We propose several research directions which we believe are novel. In the following note we go over our theoretical insights, as well as providing an analysis of possible implementation details whenever extensive enumeration is required.

## Definitions and Known Results

Given integers $L_1, L_2, ... , L_n$ we define the associated **$\bold{E_n}$ Polynomial** recursively as:

$$
E_1(x; L) := x^2 - L \\
E_n(x; L_1, ..., L_n) := E_1(E_{n-1}(x; L_1, ..., L_{n-1}); L_n) \\
$$

An **$\bold{E_n}$ Problem** is the problem of finding such integers so the polynomial $E_n(x; L_1, ..., L_n)$ factors into linear factors over the integers. We will therefore call $(L_1, L_2, ... , L_n)$ an **$\bold{E_n}$ solution**.

- [ ] Define roots.
- [ ] Define **K-th litter condition**
- [ ] Cover Bremner.
- [ ] Cover the other paper.

## Strategy - Solution Construction

If we find two $E_n$ solutions which have identical  $L_1, ... , L_{n-1}$ coefficients, then we can construct a $E_{n+1}$ solution. We can apply this insight as a strategy for finding $E_5$ solutions by first finding many $E_4$ solutions, and then find a pair of solutions which share $P$, $Q$ and $R$.

## The Connection to Gaussian Integers

Since the 1-st Litter condition for all root pairs implies:
$$
r_{2k}^2 + r_{2k+1}^2 = 2L_1 = 2P\\
$$

we can view the root pairs as a Gaussian integers, under the constraint that they are of equal norm. When searching for roots that satisfy higher order Litter conditions, we can confine our search to the comples circle with a radius $2P$. By fixing $P$ and examining it's factorization, we can:

1. Count the number of different points with norm $2P$ (namely the sum of two squares function).

- [ ] Use conventional notation for this function.

1. Efficiently ($O(n)$) enumerate over all points (see Appendix).

Furthermore, even without a concrete $2P$ at hand, we can parameterize sets of norm-like Gaussian integers using properties of the Gaussian integers (norm-multiplicity), as we will show in the next section.

## Joint Parameterization of $2^n$ norm-like Gaussian Integers  

In this section we will show that it is possible to parameterize any set of $2^n$ Gaussian integers sharing the same norm using $2^{2^n}$ variables. Such parameterization can be helpful in finding $E_{n-2}$ solutions. It should be noted that we are interested in essentially different solutions, that means we pay special attention to the common factors and order of the roots.  

Let's start by assuming that we have a set of gaussian numbers all with norm $n$ which factors as:

$$
n = 2^{a_0}p_1^{2a_1}...p_r^{2a_r}q_1^{b_1}...q_s^{b_s}
$$

- The q's denote odd rational primes which decompose in the gaussian integers ($p = 1 \mod 4$). We shall denote their decomposition as $q_i = \mathrm{q}_i\overline{\mathrm{q}}_i = (s_i + it_i)(s_i - it_i)$.
- The p's denote odd rational primes which remain primes in the gaussian integers ($p = 3 \mod 4$). They appear with even multiplicity because TODO. The existance of such factor in the norm directly implies a common factor to all roots, therefore in most cases we can assume w.l.o.g. that our norm is free of such primes.
- The factor $2$, contributes either a global scaling of all roots (if the power is even) or a simple mutation of the solutions (rotation by 45 degrees). In either cases it will not contribute an essential change to the roots so we will ignore it as well.

Therefore any gaussian integer $x = a + ib$ whose norm is $a^2 + b^2 = n$ must satisfy (up to units $\pm 1, \pm i$ which can be safely ignored):

$$
x = \prod _{i=1}^s \mathrm{q}_i^{e_i} \overline{\mathrm{q}}_i^{b_i - e_i}
$$

Therefore be characterized by a vector $\bold{e} = (e_1, ...,  e_s)$ where $0 \le e_i \le b_i$.

For the sake of simplicity, in the next sections we will treat only the case where $b_i = 1$ so we can view $\bold{e}$ as a binary vector. Is is worth noting that in that sense, taking the NOT of certain entries of that vector is equivalent of computing the conjugate of the corresponding $\mathrm{q}_i$. The generalization beyond the binary is merely a technical obstruction.

### The $n=1$ Case and Application to the $E_3$ Problem

Let $x = a + ib, y = c + id$ be two norm-like gaussian integers with norm $n$. Following the characterization described above, we can associate $x$ with a vector $e_x$ and $y$ with $e_y$. We can associate a gaussian integer $g = m + in$ with the entries where $e_x$ and $e_y$ are identical and a gaussian integer $h = p + iq$ where they differ. We can then write:

$$
x = gh = mp - nq + i(np - mq) \\
y = g\overline{h} = mp + nq + i(np + mq) \\
$$

Which results in the known parameterization of the diophantine equation $a^2 + b^2 = c^2 + d^2$. [3] Completes the discussion in the context of $E_3$ problems. It turns out however, that the 1-st Litter condition is the only constraints required for the fulfillment of the $E_3$ problem, therefore every instance of the parameterisation is a valid $E_3$ solution.  

### The $n=2$ Case and Application to the $E_4$ Problem

#### Parameterizing 4 norm-like gaussian integers

We now consider the case where where we have four norm-like gaussian integers $n_1, n_2, n_3, n_4$ whose parameterization we seek. Same as before, we associate each with a binary vector $e_1, e_2, e_3, e_4$. The task of finding a concise parameterization now translates into finding disjoint subsets of indices, s.t. for each subset has a distinct "syndrome" over each of the four vector. Each subset has a corresponds to a "base form" which is a gaussian integer, defined by the product of gaussian decomposable primes composing $e_1$ at the subset's indices. The "syndrome" of a given subset refers to for which of $e_2, e_3, e_4$ do the base for appears conjugated. There are exactly 8 different syndromes, which we will denote by the numbers 0-7. Therefore we require exactly 8 subsets to fully parameterize our for integers. We will denote the resulting gaussian integers corresponding to each subset (and syndrome) with $X_i = s_i + i t_i$. We can therefore construct the following parameterization for 4 norm-like gaussian integers:

$$
n_1 = A + iB = X_0 X_1 X_2 X_3 X_4 X_5 X_6 X_7 \\
n_2 = C + iD = X_0 X_1 X_2 X_3 \overline{X_4 X_5 X_6 X_7} \\
n_3 = E + iF = X_0 X_1 \overline{X_2 X_3} X_4 X_5 \overline{X_6 X_7}\\
n_4 = G + iH = X_0 \overline{X_1} X_2 \overline{X_3} X_4 \overline{X_5} X_6 \overline{X_7} \\
$$

- Each of the coefficients $A,B,C,D,E,F,G,H$ of the resulting $n_i$ is a degree 8 polynomial in 16 variables.

### The General Case

While we won't go into the details here, it is easy to see how we can generalize our parameterization to any number of norm-like gaussian integers (not only powers of two). Specifically we could (in principle) write down the 2-nd and 3-rd Litter conditions for the eight gaussian integers required for the satisfaction of the $E_5$ problem (it would require $2^8 = 256$ variables and the resulting constraints would be of degree no less than $2^7 \cdot 4 = 512$).

## Formulation of the $E_4$ Problem using a single constraint

Given the parameterization of One can be easily convinced that the task of satisfying the 2-nd Litter condition, the "Q condition":

$$
A^2B^2 + C^2D^2 = E^2F^2 + G^2H^2
$$

Now translates to finding the zeros of a degree 32 polynomial in 16 variables. To the best of our knowledge, this is the first reduction of the $E_4$ problem into a single constraint problem (as opposed to the original formulation which involded 4 equations).

Investigating this polynomial is a difficult task due to it's size. The actual polynomial expression is far too big to be displayed here, and not even easily manipulated in a computer program. Therefore we need ways to simplify this result in that will still be useful, hence the following observation: our quadruplet will not always contain all available syndromes $X_0, ... , X_7$, so we can discuss the different "syndrome patterns" different solutions might have, not all necessarily the complete one. For e.g. a quadruplet resulting as:

$$
n_1 = X_0 X_1 X_4 X_7 \\
n_2 = X_0 X_1 \overline{X_4 X_7} \\
n_3 = X_0 X_1 X_4 \overline{X_7}\\
n_4 = X_0 \overline{X_1} X_4 \overline{X_7} \\
$$

Will be said to have a "${\{0,1,4,7\}}$" syndrome pattern. It is therefore natural to ask two questions:

- What syndrome patterns can we identify in the $E_4$ solutions we have discovered so far (either via enumeration or via the families of [3])?
- Are there some syndrome patterns that we can rule out altogether? are some more promising search grounds then others?

One should also note that in the context of the $E_4$ problem, the integers are divided into a pair of pairs:

$$
((n_1, n_2),(n_3, n_4))
$$

And the Q condition is agnostic to order within the pairs and in between them. This imposes an equivalence on some of the syndrome pattenrs.  

# More Ideas

- Providing an estimate for the $E_{4/5}$ solution **density** can be a powerful tool in estimating the extent of enumeration needed to find them.
- Imposing generalized linear constraints on the roots and searching for solution families with techniques similar to [3].
- Parameterizing the 2-nd Litter condition as a gaussian integer.
  - Both in standard form
  - and in "centralized form"
- Using hurwitz quaternions for the 2-nd Litter condition.

# Appendix

## Characterization of Known $E_4$ solution to syndromes

- [ ] Bremner families
- [ ] Sporadic solutions

## Efficient enumeration over all Gaussian integers with norm $2P$

# References

- [1] <https://www.researchgate.net/publication/220975993_Few_Product_Gates_But_Many_Zeros>
- [2] <https://mathworld.wolfram.com/SumofSquaresFunction.html>
- [3] <https://projecteuclid.org/journalArticle/Download?urlid=em%2F1243429952>

# Archive

We are investigating chains polynomials created by iterated squaring and addition (normal d-gems according to Brochert et. al).

$$
E(x; L) := x^2 - L \\
E(x; L_1, L_2) := E(E(x; L_1); L_2) = (x^2 - L_1)^2 - L_2 \\
E(x; L_1, L_2, L_3) := E(E(E(x; L_1); L_2); L_3) = ((x^2 - L_1)^2 - L_2)^2 - L_3 \\
E(x; L_1, ..., L_n) := E(E(x; L_1, ..., L_{n-1}); L_n) \\
$$

We seek to find $E(x; L_1, ..., L_n)$ which have all of their $2^n$ roots being integers.

If $E(x; L_1, ..., L_n)$ splits over the integers, then:

1. $L_n$ is a perfect square.
2. $E(x; L_1, ..., L_{n}) = E(x; L_1, ..., L_{n-1} - \sqrt{L_n})E(x; L_1, ..., L_{n-1} + \sqrt{L_n})$. Therefore $E(x; L_1, ..., L_{n-1} - \sqrt{L_n})$ and $E(x; L_1, ..., L_{n-1} + \sqrt{L_n})$ both split over the integers.

The process continues all the way "down" until we reach term of the form:

$$
E_{\nu}(x) = E\left(x; L_1 \pm \sqrt{L_2 \pm \sqrt{... \pm \sqrt{L_n}}}\right)
$$

Where $\nu$ encodes the $2^{n-1}$ different sign permutations possible. The above terms split iff $r_{\nu}^2 := L_1 \pm \sqrt{L_2 \pm \sqrt{... \pm \sqrt{L_n}}}$ is a perfect square as well.

We can then work out a full binary tree with $2^{n-1}$ leafs denoted $r_0, ..., r_{2^{n-1} - 1}$.

$$
r_0^2 = L_1 + \sqrt{L_2 + \sqrt{... + \sqrt{L_n}}} \\
r_1^2 = L_1 - \sqrt{L_2 + \sqrt{... + \sqrt{L_n}}} \\
r_2^2 = L_1 + \sqrt{L_2 - \sqrt{... + \sqrt{L_n}}} \\
r_3^2 = L_1 - \sqrt{L_2 - \sqrt{... + \sqrt{L_n}}} \\
... \\
r_{2^{n-1}-1}^2 = L_1 - \sqrt{L_2 - \sqrt{... - \sqrt{L_n}}} \\
$$

Writing this down this way reveals the following relations:

1.

$$
r_{2k}^2 + r_{2k+1}^2 = 2L_1 \\
$$

1.

$$
(r_{4k}^2 - r_{4k+1}^2)^2 + (r_{4k + 2}^2 - r_{4k+3}^2)^2 = 8L_2\\
$$

1.

$$
((r_{8k}^2 - r_{8k+1}^2)^2 - (r_{8k + 2}^2 - r_{8k+3}^2)^2)^2 + ((r_{8k + 4}^2 - r_{8k+5}^2)^2 - (r_{8k + 6}^2 - r_{8k+7}^2)^2)^2 = 32L_3\\
$$

1. It can be worked out easily how to continue this sequence.

Toying around with the formulas for $L_1, L_2$ (losing the generic $k$ notation for brevity):

$$
(r_0^2 + r_1^2)^2 + (r_2^2 + r_3^2)^2 = r_0^4 + r_1^4 + r_2^4 + r_3^4 + 2(r_0^2r_1^2 + r_2^2r_3^2) = 8L_1^2 \\
(r_0^2 - r_1^2)^2 + (r_2^2 - r_3^2)^2 = r_0^4 + r_1^4 + r_2^4 + r_3^4 - 2(r_0^2r_1^2 + r_2^2r_3^2) = 8L_2 \\
$$

Therefore:

$$
r_0^4 + r_1^4 + r_2^4 + r_3^4 = 4(L_1^2 + L_2) = Q_a \\
r_0^2r_1^2 + r_2^2r_3^2 = 2(L_1^2 - L_2) = Q_m \\
$$

And:

$$
(r_0^2 - L_1)^2 + (r_2^2 - L_1)^2 = \frac{1}{4}(r_0^2 - (2L_1 - r_0^2))^2 + \frac{1}{4}(r_2^2 - (2L_1 - r_2^2))^2 = 2L_2 \\
$$

### The endgame (one possible approach)

A solution to $E_3$ consists of **two pairs** of integers which satisfy the "$L_1$" condition:
$$
a^2 + b^2 = c^2 + d^2 = 2L_1 \\
$$

A solution to $E_4$ consists of **two pairs of pairs** of itegers, where all four pairs of integers satisfy the same $L_1$ condition, and in addition to that, the pairs of pairs satisfy the $L_2$ condition.

Going up another notch, a solution to $E_5$ will consist of **two pairs of pairs of pairs...**. The point here is if we spot two $E_4$ solutions which share $L_1$, $L_2$ and $L_3$, we are done...

A neccesary step for this approach is to efficiently parameterize and generate $E_4$ solutions (possibly by taking Bremner's Elliptic Curves work a step forward, or by efficient numerical search). The hard part would be to find many solutions, searching for compatible "twin" solutions from the set of solutions found would hopefully be easy.

## Starting with $E_2 = E(x; P, Q)$

In order for  $(x^2 - P)^2 - Q$ to have integer roots, it must be that $Q$ is a square of some integer, so $\sqrt{Q} \in \mathbb{Z}$.

We than have $E_2(x) = (x^2 - P - \sqrt{Q})(x^2 - P + \sqrt{Q})$. In order for $E_2$ to split to four different integer roots, we must have that both $P + \sqrt{Q}$ and $P - \sqrt{Q}$ will be perfect squares.

Let's denote:

$$
a^2 = P + \sqrt{Q} \\
b^2 = P - \sqrt{Q} \\
$$

So w.l.o.g we are free to restrict ourselves to $0 \le b \le a$.

### Observations

$$
a^2 + b^2 = 2P  \\
a^2 - b^2 = 2\sqrt{Q} \\
$$

so $2P$ is a sum of two perfect squares.

### Characterization of solutions

Finding pairs of $0\le b \le a$ s.t. $a^2 + b^2 = 2P$ will yield a solution:

$$
P = (a^2 + b^2) / 2 \\
Q = (a^2 - b^2)^2 / 4 \\

r_i = \{a, -a, b, -b\} \\
$$

### How to find solutions

Actually here it is trivial. Every choice of $0 \le b \le a$ will yield a valid $P,\ Q$ pair and an $E_2(x)$ poly that splits over the integers.

## Specifying to $E_3 = E(x; P, Q, R)$

In the $n=3$ case it only suffices to find two pairs: $(a,b),\ (c,d)$ s.t.:
$$
a^2 + b^2 = c^2 + d^2 = 2P \\
$$

This alone will suffice to construct a $E_3$ polynomial with:

$$
P = \frac{a^2 + b^2}{2} \\
\ \\
Q = \frac{1}{8}((a^2 - b^2)^2 + (c^2 - d^2)^2) \\
\ \\
R = \frac{1}{64}\left((a^2 - b^2)^2 - (c^2 - d^2)^2\right)^2 \\
$$

An observation that we see here that will hold true in more complex cases is: The last two "L coefficients" do not pose constraints, but rather they are calculated by constraints imposed by the other coefficients. So here the only constraint is imposed by $P$, and the $Q, R$ coeffs are calculated.

Another observation is, we can order the the pairs as: $0 \le b \le a,\ 0 \le d \le c,\ d \le b$.

## Quadrics

The theory of quadrics (TODO: citation needed) deals with finding rational roots to degree 2 polynomial equations. An important degree 2 equation we are interested in is:

$$
a^2 + b^2 = c^2 + d^2
$$

Since a solution to it yields an $E_3$ solution.

The "quadrics trick" works as following:

1. The assignment $\bold{r_0} = (0,1,0,1)$ is a valid solution to the equation.
2. let's consider the line $\bold{r}(t) = \bold{r_0} + \bold{v}t$ where $\bold{v}$ is a 4-vector with rational coordinates. Now we plug $\bold{r}(t)$ back to the equation:

$$
(tv_a)^2 + (1+tv_b)^2 = (tv_c)^2 + (1+tv_d)^2
$$

Rearranging:

$$
t^2 (v_a^2 + v_b^2 - v_c^2 - v_d^2) + 2t(v_b - v_d) = 0
$$

It is easy to see that $t=0$ is a solution (which shouldn't surprise us since $\bold{r_0}$ is a solution), but what is the other solution?

$$
t = \frac{1}{2}\frac{v_d - v_b}{v_a^2 + v_b^2 - v_c^2 - v_d^2}
$$

Plugging it back to $\bold{r}(t)$ we get:

$$
(a, b, c, d) = \frac{1}{2}\frac{1}{v_a^2 + v_b^2 - v_c^2 - v_d^2}(v_a (v_d - v_b), 2v_a ^2 - 2v_c^2 - 2v_d^2 + v_b^2 + v_b v_d, v_c (v_d - v_b), 2v_a ^2 - 2v_c^2 - v_d^2 + 2v_b^2 - v_b v_d )
$$

So given any rational direction vector $\bold{v} = (v_a, v_b, v_c, v_d)$ we get a rational solution to the equation. And given a rational solution, we can take out the common denominator and get an integer solution.

Seeing this the other way around, any integer solution $\bold{r}$ will be reached this way since it corresponds to a rational direction vector $\bold{v} = (a, b-1, c, d-1)$.

## Gaussian integers

Recall that for $2 < n$ we have the constraint:

$$
r_{2k}^2 + r_{2k+1}^2 = 2P \\
$$

The Guassian integers $\mathbb{Z}[i]$ is an Euclidian domain with norm $N(a + bi) = a^2 + b^2$.

Therefore the constraint can be rephrased as: **Finding $2^{n-1}$ Guassian integers $r_{2k} + ir_{2k+1}$ with equal norm**. Since $\mathbb{Z}[i]$ is UFD (unique factorization domain), each Gaussian integer factorizes uniquely (up ot the units $\pm 1, \pm i$) to primes. The guassian primes are:

1. $\mathfrak{p}_{2,\pm} = 1 \pm i$
2. $p$ if $p$ is a prime over the integers and $p = 3\mod 4$.
3. $\mathfrak{p}_{p,\pm}$ if $p$ is a prime over the integers and $p = 1\mod 4$.  

We will denote $p$ odd primes with $p = 3\mod 4$ (true Gaussian primes) and $q$ odd primes with $q = 1\mod 4$ (Gaussian decomposable primes).

It is attractive to think of the constraint problem as starting with a given $P$ (which is $L_1$) and finding all Gaussian integers with norm $2P$.

Now we think of a guassian integer $a + ib$ whose norm is $2P$. We can relate the factorization of $a + ib$ to the factorization of $2P$.

$$
a + ib = P_2 P_g P_d
$$

Where:

1. $P_2$ is the part that involves $\mathfrak{p}_{2,\pm}$ and their powers.
   1. It either contains a power of $2$ or powers of  $\mathfrak{p}_{2,+}$  or $\mathfrak{p}_{2,-}$ or a mix.
   2. Since $\mathfrak{p}_{2,+}^3 = -2\mathfrak{p}_{2,-}$ (and vice versa),
2. $P_g$ is the part involving Gaussian primes.
   1. Those primes are irrelevant since they just scale both coordinates.
   2. The resulting norm ($2P$) will have them as squares. As there is no way to get an odd power of them in the norm, we can discard norms with odd powers of guassian primes.
   3. We will always assume our numbers do not contain Gaussian primes. Since if they scale  one of the guassian integers, they appear in the norm, and thus scale all other integers with the same norm.
3. $P_d$ is the part involving decomposable primes.
   1. Each decomposable prime $q$ factorizes into two conjugate numbers $\mathfrak{p}_{q,+}$ and $\mathfrak{p}_{q,-}$.
   2. We can denote the decomposition by two numbers: $e_q,j_q$ where $e_q$ is the total power of $q$ in the norm, and $-e_q \le j_q \le e_q$. ($j_q - e_q$ is constrained to be even).
   3. $P_{q, e_q, j_q} = q^{(e_q - |j_q|)/2} \mathfrak{p}_{q,sign(j_q)}^{|j_q|}$

We can therefore parameterize our search as follows, let $2P$ be some number. We factorize it over the integers.

$$
2P = \prod _p p^{e_p}
$$

It is helpful to split the factorization into Gaussian and decomposable primes (and two):

$$
2P = 2^{e_2}\prod _p p^{e_p} \prod _q q^{e_q}
$$

1. We are free to disregard the Gaussian prime part.
2. As for the "two" part, we can address only the $2$ and $2^2$ cases, since $2^n$ for $n>2$ will still come from at most one "interesting" $\mathfrak{p}_{2,+}$ and the rest are boring $2$ scalings.

$$
2P = 2^{\{1,2\}} \prod _q q^{e_q}
$$

The Gaussian integers that give rise to such norms can be parameterized as:

$$
a + ib = \{2, \mathfrak{p}_{2,+}\} \prod _q P_{q, e_q, j_q} = \prod _q q^{(e_q - |j_q|)/2} \mathfrak{p}_{q,sign(j_q)}^{|j_q|}
$$

To sum it up, a generic norm is characterized by the multiplicities of the prime factors:

$$
2P \leftrightarrow (e_2, e_5, e_{13}, ..., e_q, ...) \\
e_2 \in \{1,2\} \\
0 \le e_q
$$

And the corresponding Gaussian integer:

$$
a + ib \leftrightarrow (e_2, (e_5, j_5), (e_{13}, j_{13}), ..., (e_q, j_q), ...) \\
e_2 \in \{1,2\} \\
0 \le e_q \\
-e_q \le j_q \le e_q \ \text{even distance} \\
$$

So the total number of options for numbers with a norm $2P$ is:

$$
\#(2P) = \prod _q (e_q + 1)
$$

Since the multiplicities are given and there are $e_q + 1$ options for each $j_q$.

Note: If we are interested only in the integers  $0 \le b \le a$ we need to divide this number by 2.

Note: This is true only if there is a prime with odd multiplicity. If there is not (all are squares), then there is an edgecase.

      Probably using the terminology of "Ideals" and "Prime Ideals" would be more accurate here. 

## Fully Parameterizing $E_3$ and drawing insight from Gaussian integers

We can rephrase the "$E_3$ problem" as finding pairs of Gaussian integers with the same norm.
Let's consider two such numbers $x = a + ib$ and $y = c + id$ whose norms satisfy $N(x)= N(y)$.

since the numbers have the same norm, they must only differ by conjugations of their prime factors (not the Gaussian primes, since they are not imaginary thus conjugating them is meaningless).

      note, conjugating the 2 part is also not interesting since it will only yield an interchange in imaginary and real coordinates. 

So looking at the "decomposable primes" part of the two numbers, we can identify the $H$, the "xor" (the part that is conjugated) and $G$, the "nxor" (the part that is the same) of the two numbers. We can denote it as follows:

$$
x = GH \\
y = G\overline{H} \\
$$

How can we use this insight to create a parameterization?

We can denote:

$$
G = m + in \\
H = p + iq \\
$$

Therefore:

$$
x = (mp - nq) + i(mq + np) \\
y = (mp + nq) + i(mq - np) \\
$$

and equivalently:

$$
a = mp - nq \\
b = mq + np \\
c = mp + nq \\
d = mq - np \\
$$

Voila: the equation at the middle left section of page 2 of Bremner.

The norm, which is equal to $2P$, is $N(G)N(H) = (m^2 + n^2)(p^2 + q^2)$

### How would we generalize this to $E_4$

#### Naive (but interesting!) approach

A naive approach would be to use three partitions, $G=m+in$, $H=p+iq$, $I=s+it$ and the four points will be:

$$
a+ib = GHI \\
c+id = GH\overline{I} \\
e+if = G\overline{H}I \\
g+ih = G\overline{H}\overline{I} \\
$$

While it has some appeal to it, it is overconstrained. Not all four-tuples of points with the same norm can be represented this way.

To take this a step further, we can examine the condition for this four-tuple to be a $E_4$ solution. On top of the four points having the same norm (The $P$, or $L_1$ condition) we also need the points to satisfy the $Q$ or $L_2$ condition which is:

$$
L_2(a,b,c,d) = L_2(e,f,g,h)\\
$$
which implies
$$
(a^2 - b^2)^2 + (c^2 - d^2)^2 = (e^2 - f^2)^2 + (g^2 - h^2)^2 \\
$$

But can also take other forms (since together with the $L_1$ constraint we can do some nice manipulations)

$$
a^4 + b^4 + c^4 + d^4 = e^4 + f^4 + g^4 + h^4 \\
\ \\
a^2b^2 + c^2d^2 = e^2f^2 + g^2h^2 \\
\ \\
(P - a^2)^2 + (P - c^2)^2 = (P - e^2)^2 + (P - g^2)^2
$$

Anyway, plugging in the explicit expression of $a,b,c,d,e,f,g,h$ into the $Q$ condition we get the following constraint

$$
32mnpq(m - n)(m + n)(p - q)(p + q)(s^2 - 2st - t^2)(s^2 + 2st - t^2) = 0
$$

And

- $m,n,p,q = 0$ yields a degenerate solution.
- $m=\pm n$, $p=\pm q$ yields a degenerate solution.
- $s^2 \pm 2st - t^2$ has no solutions in the integers (try to solve modulu 3, you get $s,t = 0 \mod 3$, but we implicitly assumed they are coprime because otherwise this solution would be a scaling of another smaller solution)

#### A better approach

Something like:

$$
a+ib = X_0 X_1 X_2 X_3 X_4 X_5 X_6 X_7 \\
c+id = X_0 X_1 X_2 X_3 \overline{X_4 X_5 X_6 X_7} \\
e+if = X_0 X_1 \overline{X_2 X_3} X_4 X_5 \overline{X_6 X_7}\\
g+ih = X_0 \overline{X_1} X_2 \overline{X_3} X_4 \overline{X_5} X_6 \overline{X_7} \\
$$

I am pretty sure that this allows us to get to all possible quadruplets but it is too heavy. Each coordinate is a polynomial of degree 8 in 16 variables. sympy just even compute the $Q$ condition, since it is a degree 32 in 16 variables and nearly all terms exist.

One can try to find "special" forms of solutions (eliminate some of the $X$s). I did it in `generic_parameterizations.ipynb` and tried to play with it.

There are many symmetries here and I am pretty sure we can do better. Also in terms of degrees of reedom. We have 16 while in reality there are only 8 (7 if you count the $P$ constraint).

## Symmetries

We are only interested in solutions that satisfy $0 \le b \le a$ which is only 1/8 of the Gaussian plane.

## Random

#### More field extension

I thought maybe $\mathbb{Q}[\sqrt{i}]$ which is a degree four extension whose norm is
$$
a^4 + 2a^2b^2 + 4a^2cd - 4abc^2 + 4abd^2 + b^4 - 4b^2cd + c^4 + 2c^2d^2 + d^4
\
$$

Might enable us to do something, not sure it does...
