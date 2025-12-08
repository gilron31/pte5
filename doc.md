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


      Probably using the terminology of "Ideals" and "Prime Ideals" would be more accurate here. 

## Fully Parameterizing $E_3$ and drawing insight from gaussian integers 

We can rephrase the "$E_3$ problem" as finding pairs of gaussian integers with the same norm. 
Let's consider two such numbers $x = a + ib$ and $y = c + id$ whose norms satisfy $N(x)= N(y)$. 

since the numbers have the same norm, they must only differ by conjugations of their prime factors (not the gaussian primes, since they are not imaginary thus conjugating them is meaningless).

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

We are only interested in solutions that satisfy $0 \le b \le a$ which is only 1/8 of the gaussian plane.

## Random

#### More field extension

I thought maybe $\mathbb{Q}[\sqrt{i}]$ which is a degree four extension whose norm is 
$$
a^4 + 2a^2b^2 + 4a^2cd - 4abc^2 + 4abd^2 + b^4 - 4b^2cd + c^4 + 2c^2d^2 + d^4
\
$$

Might enable us to do something, not sure it does...