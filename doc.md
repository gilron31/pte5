# Deriving basic results

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

2.

$$
(r_{4k}^2 - r_{4k+1}^2)^2 + (r_{4k + 2}^2 - r_{4k+3}^2)^2 = 8L_2\\
$$

3.

$$
((r_{8k}^2 - r_{8k+1}^2)^2 - (r_{8k + 2}^2 - r_{8k+3}^2)^2)^2 + ((r_{8k + 4}^2 - r_{8k+5}^2)^2 - (r_{8k + 6}^2 - r_{8k+7}^2)^2)^2 = 32L_3\\
$$

4. It can be worked out easily how to continue this sequence.

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

## Gaussian integers

Recall that for $2 < n$ we have the constraint:

$$
r_{2k}^2 + r_{2k+1}^2 = 2P \\
$$

The Guassian integers $\mathbb{Z}[i]$ is an Euclidian domain with norm $N(a + bi) = a^2 + b^2$.

Therefore the constraint can be rephrased as: **Finding $2^{n-1}$ Guassian integers $r_{2k} + ir_{2k+1}$ with equal norm**. Since $\mathbb{Z}[i]$ is UFD (unique factorization domain), each gaussian integer factorizes uniquely (up ot the units $\pm 1, \pm i$) to primes. The guassian primes are:

1. $\mathfrak{p}_{2,\pm} = 1 \pm i$
2. $p$ if $p$ is a prime over the integers and $p = 3\mod 4$.
3. $\mathfrak{p}_{p,\pm}$ if $p$ is a prime over the integers and $p = 1\mod 4$.  

We will denote $p$ odd primes with $p = 3\mod 4$ (true gaussian primes) and $q$ odd primes with $q = 1\mod 4$ (gaussian decomposable primes).

It is attractive to think of the constraint problem as starting with a given $P$ (which is $L_1$) and finding all gaussian integers with norm $2P$.

Now we think of a guassian integer $a + ib$ whose norm is $2P$. We can relate the factorization of $a + ib$ to the factorization of $2P$.

$$
a + ib = P_2 P_g P_d
$$

Where:

1. $P_2$ is the part that involves $\mathfrak{p}_{2,\pm}$ and their powers.
   1. It either contains a power of $2$ or powers of  $\mathfrak{p}_{2,+}$  or $\mathfrak{p}_{2,-}$ or a mix.
   2. Since $\mathfrak{p}_{2,+}^3 = -2\mathfrak{p}_{2,-}$ (and vice versa),
2. $P_g$ is the part involving gaussian primes.
   1. Those primes are irrelevant since they just scale both coordinates.
   2. The resulting norm ($2P$) will have them as squares. As there is no way to get an odd power of them in the norm, we can discard norms with odd powers of guassian primes.
   3. We will always assume our numbers do not contain gaussian primes. Since if they scale  one of the guassian integers, they appear in the norm, and thus scale all other integers with the same norm.
3. $P_d$ is the part involving decomposable primes.
   1. Each decomposable prime $q$ factorizes into two conjugate numbers $\mathfrak{p}_{q,+}$ and $\mathfrak{p}_{q,-}$.
   2. We can denote the decomposition by two numbers: $e_q,j_q$ where $e_q$ is the total power of $q$ in the norm, and $-e_q \le j_q \le e_q$. ($j_q - e_q$ is constrained to be even).
   3. $P_{q, e_q, j_q} = q^{(e_q - |j_q|)/2} \mathfrak{p}_{q,sign(j_q)}^{|j_q|}$

We can therefore parameterize our search as follows, let $2P$ be some number. We factorize it over the integers.

$$
2P = \prod _p p^{e_p}
$$

It is helpful to split the factorization into gaussian and decomposable primes (and two):

$$
2P = 2^{e_2}\prod _p p^{e_p} \prod _q q^{e_q}
$$

1. We are free to disregard the gaussian prime part.
2. As for the "two" part, we can address only the $2$ and $2^2$ cases, since $2^n$ for $n>2$ will still come from at most one "interesting" $\mathfrak{p}_{2,+}$ and the rest are boring $2$ scalings.

$$
2P = 2^{\{1,2\}} \prod _q q^{e_q}
$$

The gaussian integers that give rise to such norms can be parameterized as:

$$
a + ib = \{2, \mathfrak{p}_{2,+}\} \prod _q P_{q, e_q, j_q} = \prod _q q^{(e_q - |j_q|)/2} \mathfrak{p}_{q,sign(j_q)}^{|j_q|}
$$

To sum it up, a generic norm is characterized by the multiplicities of the prime factors:

$$
2P \leftrightarrow (e_2, e_5, e_{13}, ..., e_q, ...) \\
e_2 \in \{1,2\} \\
0 \le e_q
$$

And the corresponding gaussian integer:

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

## Symmetries

We are only interested in solutions that satisfy $0 \le b \le a$ which is only 1/8 of the gaussian plane.

# Unordered ideas

1. Multiplying by $\mathfrak{p}_{2,+} = 1 + i$ is like transforming $(a, b) \rightarrow (a + b, a - b)$
2. Can we prove that the existence of a double root makes it impossible to find solutions?

# Archive

$$
P_{00} = G_0 \overline{G_1} G_2 \overline{G_3} \\
P_{01} = G_0 \overline{G_1} \overline{G_2} G_3 \\
$$

And the second by

$$
P_{10} = G_1 G_1^c\\
P_{11} = G_1 \overline{G_1^c}\\
$$

## Draft 1

The first pair is parameterized by
$$
P_{00} = G_0 G_0^c\\
P_{01} = G_0 \overline{G_0^c}\\
$$

And the second by

$$
P_{10} = G_1 G_1^c\\
P_{11} = G_1 \overline{G_1^c}\\
$$

Let's denote

$$
G_{i,j} = G_i \cap G_j \\
G_{i,j^c} = G_i \cap G_j^c \\
$$

There are thus four components:

$$
G_{0,1} \\
G_{0,1^c} \\
G_{0^c,1} \\
G_{0^c,1^c} \\
$$

$$
G_0 = G_{0,1} \cup G_{0,1^c} \\
G_1 = G_{0,1} \cup G_{0^c,1} \\
G_0^c = G_{0^c,1} \cup G_{0^c,1^c} \\
G_1^c = G_{0,1^c} \cup G_{0^c,1^c} \\
$$

Let's rewrite the pairs:

$$
P_{00} = G_{0,1} G_{0,1^c} G_{0^c,1} G_{0^c,1^c}\\
P_{01} = G_{0,1} G_{0,1^c} \overline{G_{0^c,1}} \overline{G_{0^c,1^c}}\\
P_{10} = G_{0,1} G_{0,1^c} G_{0^c,1} G_{0^c,1^c}\\
P_{11} = G_{0,1} G_{0,1^c} G_{0^c,1} G_{0^c,1^c}\\

$$
