# Sources — Key Insights

Synthesis of the four source PDFs, written while restructuring `paper.md`.
Purpose: (a) resolve the open literature questions in `TODO.md` (W0), (b) fix the
correspondence/novelty issues from `critical_review.md`, (c) record the exact
numbers that go into the paper, and (d) index further references worth obtaining.

All numeric claims below were re-derived in `drafts/verify_claims.py` (run it to
reproduce). **Headline correction: the review's E₅ bound dropped a prime — see §4.**

---

## 0. Notation dictionary (each paper uses different conventions)

| Object | This paper | Bremner [Bre08] | BMR [BMR09] | CSMV [CSMV23] | TLT [TLT26] |
|---|---|---|---|---|---|
| The tower | $E_n$, degree $2^n$ | $(((X^2-P)^2-Q)^2-R)^2-S^2$ ($=E_4$) | skew / normalized $2^n$-gem | — | — |
| Tower coeffs | $L_1,L_2,L_3,L_4,\dots$ | $P,Q,R,S^2$ $=L_1,L_2,L_3,L_4$ ($S=\sqrt{L_4}$) | $\gamma_i=-L_i$ | — | — |
| Tree/level conditions | $L_n$ conditions | (equal power sums) | **Litter conditions** (Def 5.2) | — | — |
| PTE solution | size $2^{n-1}$, degree $2^{n-1}{-}1$ | 8 pairs $\pm A..\pm H$ | size $n$, degree $k$; ideal $k{=}n{-}1$ | size $n$, degree $m{=}n{-}1$ | degree $k{=}3$, $n{=}4$ |
| Recursively-symmetric subclass | — | — | **sym-perfect** (Def 3.6) | — | — |
| $L_1$ condition | $r_{2k}^2{+}r_{2k+1}^2{=}2L_1$ | $A^2{+}B^2{=}\dots$ | first litter condition | — | $x^2{+}y^2{=}u^2{+}v^2$ |

**Size ↔ degree.** An $E_n$ solution yields a symmetric ideal PTE solution of
**size $2^{n-1}$** (= number of elements per side) and **degree $2^{n-1}-1$**.
So $E_4\to$ size 8 / degree 7; $E_5\to$ size 16 / degree 15. The proposal's
"degree $n=15$" is CSMV's "size 16." Use **size** to avoid confusion.

---

## 1. Bremner 2008 — this is reference [3], and it is the $E_4$ paper

**A. Bremner, "When Can $(((X^2-P)^2-Q)^2-R)^2-S^2$ Split into Linear Factors?",
*Experimental Mathematics* 17:4 (2008), 385–390.** (Project Euclid `em/1243429952`.)

- **Identity of the object.** Bremner's polynomial *is* our $E_4$, with
  $(P,Q,R,S)=(L_1,L_2,L_3,\sqrt{L_4})$. He splits it into 16 linear factors
  $(X^2-A^2)\cdots(X^2-H^2)$. So the "first family / second family / sporadic
  Table 1" that `paper.md` analyzes are **all Bremner's**.
- **He poses $E_5$ as open**: "Whether there can exist examples of such
  identities for five or more nested squares is an open question." Our paper is
  the direct sequel to this sentence.
- **The variety $V$** (his §2): the $E_4$ locus is a **degree-32 threefold in
  $\mathbb P^7$**, cut by $A^2{+}B^2{=}C^2{+}D^2{=}E^2{+}F^2{=}G^2{+}H^2$ (three
  quadrics, the $L_1$ conditions) and $A^2B^2{+}C^2D^2=E^2F^2{+}G^2H^2$ (one
  quartic, the $L_2$/$L_{2,b}$ condition). Equivalently, equal power sums of
  degree 2, 4, 6 for $\{A,B,C,D\}$ vs $\{E,F,G,H\}$. **This is the same object as
  our §5 $W$-reduction** → the "first single-constraint reduction" claim MUST be
  hedged: Bremner already has the variety as a complete intersection; our
  contribution is a *single polynomial* inside an explicit *multiplicative
  (Gaussian) parameterization* that trivializes the three quadrics.
  - Degree $2\cdot2\cdot2\cdot4=32$, dimension $7-4=3$ ⇒ confirms the §8 geometry
    claim "type $(2,2,2,4)$ complete intersection in $\mathbb P^7$."
  - Curious identity he notes: $(A^4-B^4)^2+(C^4-D^4)^2=(E^4-F^4)^2+(G^4-H^4)^2$
    — this is the $L_{2,d}$/fourth-power structure of our §5.
- **Both families satisfy the linear condition $A-C=E-G$** (his (2-3)); the second
  starred solution also $A{+}B{+}C{+}D=E{+}F{+}G{+}H$ (his (2-4)).
  - **First family**: intersect $V$ with $A-C=E-G$ → elliptic curve
    $Y^2=(X-4)(X-3)(X+6)$, rank 1, generator $P=(2,4)$; explicit degree-(2,3)
    parameterization of $A,\dots,H$ in $(x,y,z)$ on the quadric
    $\mu^2+\phi^2=2\nu^2+10\lambda^2$.
  - **Second family**: genus-1 curve $\mathcal C: 3x^2y-4xy^2+y^3+(x^2-14xy+9y^2)z+10(x+y)z^2=0$,
    birational to $Y^2=X^3-X^2-8X+112$, rank 1, generator $(12,40)$.
  - **These symbolic forms unblock R6** (family intersection): solve
    $L_1(k)=L_1'(m),L_2=L_2',L_3=L_3'$ over the two curves.
- **Solutions come in pairs** (his p.4): $(A,B,C,D,E,F,G,H)\mapsto(A{+}B,A{-}B,C{+}D,C{-}D,E{+}F,E{-}F,G{+}H,G{-}H)$.
- **Origin**: Dilcher [Di00] parameterized $k=3$ (our $E_3$) and *claimed $k=4$
  impossible*; Symes/Crandall–Pomerance found a $k=4$ example, refuting it;
  Bremner then produced the two infinite families.

## 2. Borchert–McKenzie–Reinhardt 2009 — resolves TWO novelty questions

**B. Borchert, P. McKenzie, K. Reinhardt, "Few Product Gates but Many Zeros",
MFCS 2009 / *Chicago J. Theoret. Comput. Sci.* (2009).** DOI 10.1007/978-3-642-03816-7_15.

- **d-gem / skew / normalized (Defs 1.1, §3).** A *skew (≡ normalized) $2^n$-gem*
  is exactly our $E_n$ tower; a $2^n$-gem is a degree-$2^n$ polynomial with $2^n$
  distinct integer roots computed with $n$ product gates. **An $E_5$ solution is a
  skew 32-gem.**
- **Complexity motivation (§1, §6).** Skew $2^n$-gems for infinitely many $n$
  would refute the **$L$-conjecture** (Bürgisser [Bu01]) and hence the
  **$\tau$-conjecture** (Blum–Cucker–Shub–Smale), which implies
  $\mathrm P_{\mathbb C}\ne\mathrm{NP}_{\mathbb C}$ (Smale's list). Factoring link
  via Lipton. This is the algebraic-complexity hook for §1.
- **Litter conditions = our tree/$L_n$ conditions (Def 5.2).** ETYMOLOGY OF
  "Litter" RESOLVED — it is BMR's term. Example 5.3 gives, for squared roots
  $a^2,\dots,h^2$: $a^2{+}b^2{=}c^2{+}d^2{=}e^2{+}f^2{=}g^2{+}h^2$ **and**
  $a^2b^2{+}c^2d^2{=}e^2f^2{+}g^2h^2$ — verbatim our $L_1$ and $L_{2,b}$.
- **The corrected correspondence (this fixes the proposal's overclaim):**
  - *Cor 3.5*: every skew $2^n$-gem gives an ideal PTE solution of size $2^{n-1}$
    (the top split into the two $E_{n-1}$ factors).
  - *Thm 3.7*: a size-$2^{n+1}$ set is the zero set of a normalized $2^{n+1}$-gem
    **iff** it is a **sym-perfect** pair (Def 3.6 — a *recursive* condition).
  - *Cor 3.8*: sym-perfect ⟹ ideal symmetric PTE, but the inclusion is **proper**.
  - ⇒ **"$E_n$ = recursively-symmetric subclass" is precisely BMR's sym-perfect.**
    The critical review's phrase is BMR's theorem; cite it, don't re-introduce it.
- **The Tarry counterexample is in BMR (p.7).** They take the symmetric ideal
  degree-7 solution $\{2,16,21,25\},\{5,14,23,24\}$ (Borwein–Ingalls; = centered
  Tarry 1913) and note $2^2{+}25^2\ne16^2{+}21^2$, so it yields only a **non-skew**
  16-gem. That inequality is exactly the failure of our $L_1$ condition ⇒ this
  solution is **not an $E_4$**. Use directly in §2.
- **BMR tried and failed** to find a 32-gem, or even a 16-gem with 16 distinct
  squares as roots (p.10). Independent corroboration that $E_5$ is hard.
- Misc: *Prop 3.1* gem scaling (homogeneity $a_i\mapsto ta_i$); *Cor 4.5* a skew
  $2^n$-gem needs exactly $n$ additive gates.
- **Novelty impact on §3.** BMR's sym-perfect (Def 3.6) already encodes the
  recursive/level-wise structure. Our tower recursion $u_k=(r_{2k}^2-r_{2k+1}^2)/2$
  must be framed as a **congruum-lift reading of BMR**; the genuinely new packaging
  is the *bottom-lift* algorithm and the meet-in-the-middle cut points, not "we
  introduce the recursion."

## 3. Coppersmith–Mossinghoff–Scheinerman–VanderKam 2023 — the arithmetic engine

**D. Coppersmith, M. Mossinghoff, D. Scheinerman, J. VanderKam, "Ideal Solutions
in the Prouhet–Tarry–Escott Problem", arXiv:2304.11254 (2023).** Dedicated to
Peter Borwein.

- **The constant.** $C_n(A,B)=\prod(x-a_i)-\prod(x-b_i)$ (degree 0 for ideal
  solutions). $C_n=\gcd$ over ideal size-$n$ solutions; $C'_n=\gcd$ over
  *symmetric* ideal size-$n$ solutions. Kleiman: $(n-1)!\mid C_n$. Rees–Smyth
  Prop 2.1 + Multiplicity Lemma 2.1 give the required-divisor machinery;
  Filaseta–Markovich give 2-adic ($2^6\mid C_8$, $2^9\Vert C_9$).
- **Table 1 (required divisors)** — the input to our ladder. Key rows:
  - $C_4=2^2\!\cdot\!3^2$, $C_8=2^6\!\cdot\!3^3\!\cdot\!5^2\!\cdot\!7^2\!\cdot\!11\!\cdot\!13$ (no extra for $C'$).
  - **$C_{16}\supseteq 2^{11}\!\cdot\!3^6\!\cdot\!5^4\!\cdot\!7^3\!\cdot\!11^2\!\cdot\!13^2\!\cdot\!17\!\cdot\!19\!\cdot\!23$; additional for $C'_{16}$: $29\!\cdot\!37\!\cdot\!41\!\cdot\!43\!\cdot\!53$.**
- **Table 2 (exhaustive symmetric searches, all NULL above size 12):**

  | size $n$ | 9 | 10 | 11 | 12 | 13 | 14 | 15 | **16** |
  |---|---|---|---|---|---|---|---|---|
  | height $H$ | 7000 | 2500 | 3500 | 1511 | 2000 | 900 | 1100 | **850** |
  | # solutions | 2 | 2 | **0** | 2 | **0** | **0** | **0** | **0** |

  This is the prior art any $E_5$ must clear: **no symmetric ideal size-16
  solution of height $\le 850$.** (Also nothing at size 11, 13, 14, 15.)
- **State of the art.** Ideal PTE solutions over $\mathbb Z$ known for sizes
  $n\le10$ and $n=12$; **open at $n=11$ and $n\ge13$.** Largest known = size 12
  (degree 11).
- **§3 local obstructions (Chebyshev/Dickson).** For a prime $p>n$ with
  $p\equiv\pm1\pmod n$, a symmetric ideal solution exists over $\mathbb F_p$, so
  $p$ **cannot** be a required divisor of $C_n$. This *partially forecloses*
  direction 5 (sharper required divisors) for such primes — but says nothing
  about primes forced by the *towered* structure specifically. Read their §3
  before claiming any new required divisor.
- **§§5–6 rings.** New solutions only at size $\le12$ in every ring: $\mathbb Z[i]$
  at sizes 10 (five, eqns 24–28) and 12 (one, eqn 29, with $i$-rotation symmetry
  $A=\bigcup_k i^k\{3{+}10i,11{+}6i,8{+}10i\}$); Eisenstein at 9 and 12 (sixfold
  $\zeta_6$ symmetry, eqn 32); $\mathbb Z[i\sqrt2]$ at 9. **Nothing above 12 in any
  ring** — bounds the ring-pivot expectation, and supplies the rotation ansatz to
  transplant to size 16.

### The divisibility ladder + corrected $E_5$ bound (verified)

Mechanism (verified on the exemplar, §"verified" below): for an $E_5$, the top
split gives $S,T$ = roots of $E_4(x)\mp\sqrt{L_5}$, so
$C_{16}(S,T)=E_4^-(x)-E_4^+(x)=-2\sqrt{L_5}$ **exactly**. Since $(S,T)$ is a
symmetric ideal size-16 solution, $C'_{16}\mid 2\sqrt{L_5}$. Therefore:

$$\boxed{\ \sqrt{L_5}\ \ge\ C'_{16}/2\ =\ 2\,437\,428\,918\,743\,498\,865\,144\,960\,000\ \approx\ 2.44\times10^{27}\ }$$

hence $L_5\ge(C'_{16}/2)^2\approx5.94\times10^{54}$. **The review's
$1.43\times10^{26}$ is wrong: it dropped the prime 17** (that value equals
$C'_{16}/(2\cdot17)$). The height bound is *separate*: an $E_5$'s $(S,T)$ has
height $=r_{\max}$, which by CSMV Table 2 must exceed 850, and since
$r_{\max}^2\in[L_1,2L_1]$ this gives $L_1>850^2/2=361250$. The review's
"root $>3400$, $L_1>5.8\times10^6$" could **not** be reproduced — the ladder alone
gives only $L_1\ge(C'_{16}/2)^{1/8}\approx2650$ (weaker than CSMV's $361250$), so
CSMV's search, not the ladder, supplies the height bound. Cascade at lower levels
(verified): $C'_8/2=151351200\mid\sqrt{L_4}$; $C'_4/2=18\mid\sqrt{L_3\pm\sqrt{L_4}}$.

## 4. Tsai–Lee–Takahashi 2026 — our $E_3$ case, done rigorously

**Y.-D. Tsai, J. Lee, F. Takahashi, "Arithmetic Symmetry in Ideal
Prouhet–Tarry–Escott Solutions", arXiv:2606.07735 (June 2026).** Motivated by
anomaly cancellation in chiral gauge theory (Lee–Takahashi–Tsai charge spectra).

- **Their symmetric degree-3 PTE = our $E_3$.** After centering ($c\in\frac12\mathbb Z$)
  and doubling, a symmetric ideal degree-3 solution is $\{\pm x,\pm y\},\{\pm u,\pm v\}$,
  which is a solution **iff $x^2+y^2=u^2+v^2$** (their Prop 2.3) — literally our
  $L_1$ condition / equal-norm Gaussian integers.
- **Counting law**: $C_{\mathrm{sym}}(H)=\frac{2\log2}{\pi^2}H^2\log H+O(H^2)$
  (centered), $N_{\mathrm{sym}}(H)=\frac{4\log2}{3\pi^2}H^3\log H+O(H^3)$ (summed
  over centers). The **log-enhancement** is the second moment of the
  sum-of-two-squares function, $\sum_{n\le X}r_2(n)^2\asymp X\log X$ (Borwein–Choi
  asymptotic with explicit constant).
- Explicit family: $(pr-1)^2+(p+r)^2=(pr+1)^2+(p-r)^2$ (from $\alpha\beta=\gamma\delta$,
  i.e. Brahmagupta) giving $\gg H\log H$ solutions.
- §6 density conjecture for the **full** degree-3 space: $N(H)=H^{3+o(1)}$;
  Appendix A relates it to the Vinogradov mean value $J_{4,3}$, raw bound
  $N(H)\ll H^{4+\varepsilon}$.
- **Impact.** Our §8 density section **must start at $E_4$ and cite TLT for
  $E_3$** — they published the $E_3$ density one month ago via the *same*
  $\mathbb Z[i]$/sum-of-two-squares mechanism. Our framework is the $E_4$/$E_5$
  generalization; the second-moment method is the tool to lift to $E_4$.

---

## 4b. Dilcher 2000 — the $E_3$ and $E_4$ reductions predate us; the famous wrong claim

**K. Dilcher, "Nested squares and evaluations of integer products", *Exp. Math.*
9:3 (2000), 369–372.** (`sources/Dil00.pdf`.) Origin of the whole nested-squares
line; motivated by Crandall's fast integer-product identity
$((x^2-85)^2-4176)^2-2880^2=(x^2-1)(x^2-49)(x^2-121)(x^2-169)$ (an $E_3$).

- **Prop 1 ($E_3$ characterization):** if $n=a_1^2+b_1^2=a_2^2+b_2^2$ then
  $((x^2-\tfrac n2)^2-\cdots)^2-\cdots=(x^2-a_1^2)(x^2-b_1^2)(x^2-a_2^2)(x^2-b_2^2)$.
  So **$E_3$ ⟺ an even $n$ with two representations as a sum of two squares** —
  our $L_1$ condition, and TLT's $x^2+y^2=u^2+v^2$, stated by **Dilcher in 2000**,
  26 years before TLT. (TLT's *density* is new; the reduction is not.)
- **§3 ($E_4$ as a single constraint):** for four representations
  $n=a_i^2+b_i^2$, the nested form of length 4 (our $E_4$) exists **iff**
  $a_1^2b_1^2+a_2^2b_2^2=a_3^2b_3^2+a_4^2b_4^2$ (his (3-2)) — this is exactly our
  $L_2$/$L_{2,b}$ condition and **is a single scalar constraint on four equal-norm
  Gaussian integers**. ⇒ **The "single-constraint reduction of $E_4$" is
  essentially Dilcher (3-2), 2000.** Retract the "first" claim; our genuine
  contribution is the *multiplicative (syndrome) parameterization* turning (3-2)
  into a polynomial identity $W(\mathbf X)=0$, plus the enumeration algorithm.
- **Prop 2 (the famous error):** Dilcher *proved* (3-2) has no solution, i.e.
  claimed **$E_4$ is impossible**. His argument fixes one factorization pattern
  ($n=n_1n_2$) and one of the three sign patterns of
  $d=a_1^2b_1^2+a_2^2b_2^2-a_3^2b_3^2-a_4^2b_4^2$; it does not cover all
  quadruples. Symes/Crandall–Pomerance produced a counterexample and **Bremner
  [Bre08] the families** — refuting it. Excellent history for §1 (a cautionary
  tale about "impossibility" claims that directly motivates our search).
- **Remark 3 (divisibility):** the $L_2$-difference $d$ is **always divisible by
  $1152=2^7\cdot3^2$**, and by $28800=2^7\cdot3^2\cdot5^2$ when $5\nmid n$. A
  concrete forced-divisor fact feeding the §7 smoothness discussion and the ladder.
- Remark 2: a product of $2^{k+2}$ linear factors $=$ nested squares of length
  $k{+}2$ plus an error term of degree $2^{k+2}-12$.

## 4c. Raghavendran–Narayanan 2019 — PTE survey (background + reference spine)

**S. Raghavendran, V. Narayanan, "The Prouhet–Tarry–Escott Problem: A Review",
*Mathematics* 7(3):227 (2019).** (`sources/RN19.pdf`, open access.) Not a research
input but the cleanest one-stop history + bibliography for §1.

- History spine for §1: Euler/Goldbach (~1750) → Prouhet 1851 → Tarry/Escott
  (1910s) → Frolov 1889 (affine invariance $Mx_i+K$) → Bastien 1913 (nontrivial
  ⇒ size $\ge$ degree $+1$) → Dorwart–Brown 1937, Chernick 1937 (ideal
  constructions $k=3..7$) → Wright 1959 (Prouhet). Prouhet–Thue–Morse construction.
- State of the art confirmed: ideal solutions known for sizes $\le10$ and $12$;
  **size 11 and $>12$ open**. Figure 1 lists the smallest known ideal solutions
  sizes 2–12 (their notation $K_n$ = our/CSMV's constant $C_n$).
- **Open-problems list (their §5) aligns with our program**: (a) does an ideal
  solution exist beyond size 12/at 11?; (b) solve mod $p^n$ for primes $\ne$ size;
  (c) **"find a lower bound on $N(k)$ that would rule out ideal solutions"** —
  exactly the spirit of our §4 lower bounds; (d) Filaseta–Markovich 2-adic
  $v_2(\bar C_9)=9$.
- Useful adjacent references it collects (full details in its bibliography):
  Choudhry 2000/2003/2017 (degree-4/5 ideal + AP solutions), Alpers–Tijdeman 2007
  (2-D PTE; **introduced classical PTE over Gaussian integers**), Prugsapitak
  2009–2013 (PTE over quadratic fields / $\mathbb Z[i]$), Gandikota–Ghazi–Grigorescu
  2016 (**NP-hardness of Reed–Solomon decoding ↔ PTE**, a complexity hook beyond
  BMR), Bremner 1981 (equal sums of fifth powers via a K3 surface).

## 5. Novelty ledger (what is ours vs. prior) — for the "we introduce" sentences

| Ingredient | Status | Attribution |
|---|---|---|
| $E_n$ ansatz / tower of nested squares | **prior** | BMR skew gems; Bremner nested squares |
| Litter / tree $L_n$ conditions | **prior** | BMR Def 5.2 |
| $E_n$ ⟺ recursively-symmetric ⊊ symmetric ideal | **prior** | BMR Cor 3.5 / Thm 3.7 / Cor 3.8 |
| Tarry/$L_1$ counterexample | **prior** | BMR p.7 (Borwein–Ingalls solution) |
| $E_4$ = degree-32 threefold (complete intersection) | **prior** | Bremner §2 (variety $V$) |
| $E_4$ ⟺ single $L_2$ constraint on 4 equal-norm reps | **prior** | **Dilcher 2000 (3-2)** — retract "first single-constraint" |
| Multiplicative (syndrome) parameterization ⇒ $W(\mathbf X)=0$ polynomial | **ours** (angle) | this work; makes Dilcher's (3-2) a polynomial identity |
| $L_{2,a,b,c,d}$ augmentations; $L_{2,d}$ fourth-power form | **ours** (minor) | this work |
| $\mathbb Z[\sqrt2]$ theorem ($L_{2,d}$ primes $\equiv\pm1\bmod8$) | **ours** (verified) | this work |
| Gaussian *multiplicative* parameterization + single $W$ | **ours** (angle) | this work; relate to Bremner's $V$ |
| Tower recursion $u_k$ / bottom-lift / MITM | **congruum-lift reading of BMR** + new algorithms | this work reads [BMR] |
| Divisibility ladder + explicit $E_5$ lower bounds | **corollary of CSMV** (verified) | imports [CSMV] Table 1 |
| $E_4$/$E_5$ density (multidegree, singular series) | **ours** ($E_3$ is TLT) | this work; cite [TLT] |
| Per-norm enumeration algorithm at scale | **ours** (needs non-squarefree fix) | this work |

Two "we introduce" claims to **retract/hedge**: the recursion (→ BMR sym-perfect)
and the single-constraint reduction (→ Bremner's variety $V$).

## 5b. Chapter 5 proved: surjectivity + single-constraint reduction

Both §5 load-bearing claims are now proven in `paper.md` (Prop 5.1, Prop 5.2) and
checked in `drafts/verify_claims.py` §6.

- **Single constraint (Prop 5.2).** Under the multiplicative parameterization
  $\Phi(X_0,\dots,X_7)$, equal norms (the $L_1$ condition) hold identically and the
  quadruplet is an $E_4$ solution **iff** $W = A^2B^2 + C^2D^2 - E^2F^2 - G^2H^2 = 0$.
  Proof: each pair is an $E_3$ (Brahmagupta/Dilcher) with second coefficient
  $L_2^{(i)} = L_1^2 - \tfrac12(\cdot)$; two $E_3$'s already sharing $L_1$ glue into an
  $E_4$ iff they also share $L_2$ (the §2 Key Property), which is exactly $W=0$. This is
  the *exhaustive* form of Dilcher (3-2); attribution unchanged.
- **Surjectivity (Prop 5.1) — for ALL admissible norms.** The review's trichotomy
  (`critical_review.md` §1) is the right frame: **S1** ("every equal-norm quadruple has
  a multiplicative description") is conceded as trivial; what we add is the
  *constructive, multiplicity-complete* version. Key idea (**token argument**): for a
  split prime $p = \mathrm{q}\overline{\mathrm{q}}$ with $v_p = b$, distribute the $b$
  prime factors as independent tokens over the 8 buckets $\times$ the
  $\mathrm{q}/\overline{\mathrm{q}}$ choice; a single token realizes **any** vector in
  $\{0,1\}^4$ of exponent-increments $(c_1,c_2,c_3,c_4)$, so $b$ tokens realize **any**
  $(c_1,\dots,c_4) \in \{0,\dots,b\}^4$ — every exponent quadruple, including the
  primitive $p^2$ case (one pair uses the real factor $\mathrm{q}\overline{\mathrm{q}}=p$,
  $c=1$; another uses $\mathrm{q}^2$, $c=2$). Inert primes and the ramified
  $2 = -i(1+i)^2$ are global common factors ($\to X_0$). CRT over primes $\Rightarrow$
  surjectivity up to units.
- **What this does NOT close.** Surjectivity is expressive completeness, not uniqueness
  or code-coverage:
  - **Uniqueness / canonization** (still a TODO): for $b \ge 2$ the token distribution
    is non-unique, so the "syndrome pattern" is a *multiset* over buckets; the
    canonization lemma (well-defined §7 census) is still owed.
  - **Enumerator coverage "S3"** (still a §6 code TODO): the *implemented* Stage-1–3
    enumerator must iterate exponent vectors $0 \le e_i \le b_i$, exactly as Prop 5.1's
    proof prescribes, or "exhaustive to $B$" keeps its asterisk. The math half is now
    discharged; the code half is not.
- **Bit convention (reproducibility).** Bucket $t = 4\beta_2 + 2\beta_3 + \beta_4$;
  $n_2$ conjugates $X_{4,5,6,7}$, $n_3$ conjugates $X_{2,3,6,7}$, $n_4$ conjugates
  $X_{1,3,5,7}$. Reconstruction verified: exemplar (squarefree) $\to$ occupied buckets
  $\{0,1,4,7\}$; primitive $5^2\cdot13\cdot17$ quadruple $\to \Phi$ reproduces inputs up
  to units. The occupied-bucket set is reference/labeling dependent — hence the need for
  canonization; it is not literally §7's $\{0,1,3,5\}$, which is computed *after* the §7
  canonization steps.

**Novelty-ledger delta** (supersedes the two `angle` rows in §5): the multiplicative
parameterization *with proved surjectivity for all admissible norms* is **ours
(proved + verified)** — constructive token argument, makes Dilcher (3-2) exhaustive; the
$W=0$ reduction inside $\Phi$ (Prop 5.2) is **ours as packaging** (Dilcher (3-2) is the
constraint, ours is the $\Phi$-identity form plus the exhaustiveness that upgrades
"necessary" to "iff").

## 5c. Critical-pass findings (correctness / scope / novelty)

A reviewer/advisor pass over the whole paper. Fixes applied unless noted.

**Correctness bug (found + fixed + guarded).**
- **§2 level-3 tree constant was wrong: `32L_3` → `128L_3`.** Verified on the
  exemplar (`LHS/L_3 = 128`, remainder 0; `verify_claims.py` §7 now asserts this).
  The level-$k$ constant is $2^{2^k-1}$ (values $2, 8, 128, \dots$), from nesting
  the §3 recursion — *not* the arithmetic $2, 8, 32$ the old text assumed. The
  recursion (§3) and the recovered-tower checks were unaffected (they use the $8L_2$
  form, which was correct); only the displayed $L_3$ line was wrong. Worth a
  self-note: **the pattern doubles the exponent each level**, so an $E_5$'s $L_4$
  condition would carry $2^{15}=32768$.

**Scope / statement precision (fixed).**
- **Abstract correspondence.** "in pairs sharing $L_1..L_{n-1}$ … size $2^{n-1}$"
  conflated two facts and mis-indexed. A *single* $E_n$ → symmetric ideal size
  $2^{n-1}$ (BMR Cor 3.5, via its top split); "pairs sharing $L_1..L_{n-1}$" is the
  Key Property (→ $E_{n+1}$, size $2^n$). Fixed to "each yields, via its top split,
  … size $2^{n-1}$", matching the (correct) §2 body.
- **§1 complexity implication.** "(hence $P_\mathbb{C}\ne NP_\mathbb{C}$)" was
  backwards: refuting the $\tau$-conjecture does **not** establish the separation —
  the $\tau$-conjecture *implies* $P_\mathbb{C}\ne NP_\mathbb{C}$, and gems *refute*
  $\tau$. Reworded so the arrow is unambiguous. **Still to verify against BMR §6 /
  Bürgisser**: the exact implication direction between the $L$-conjecture and the
  $\tau$-conjecture (paper now just says gems refute both; safe).
- **§4 required-divisor definition.** "gcd over all size-16 solutions" is ill-posed
  when none is known. Fixed to "largest integer *proven* to divide $C_n(A,B)$ for
  every solution — a required divisor, established by local/2-adic arguments." The
  bound only needs the required-divisor property. Also trimmed the "should be
  re-checked before publication" internal note in §4 to a clean remark.

**Progress this pass (chapters 5, 6, 8).**
- **§5 Proposition 5.4 (degenerate locus of $W$)** — done: genuine $E_4$ ⟺ $W=0$
  **and** the 8 coordinates nonzero & pairwise-distinct ($\Delta\ne0$); $D=\{W=0\}\cap\{\Delta=0\}$
  is proper closed (contains $X_t=0$, zero-coordinate, and coincidence loci).
  Discharges the "degenerate locus" half of the §5 hygiene TODO; **canonization
  lemma still open** (the multiset-syndrome case for $b\ge2$).
- **§6 rewritten** into a full search design: 6.1 per-norm enumeration (with the
  non-squarefree Stage-2 branching = Prop 5.1), 6.2 two regimes + **Prop 6.1
  (exclusion)**, 6.3 top-collision, 6.4 bottom-lift, 6.5 ladder pre-sieve (with the
  completeness-is-a-theorem caveat), 6.6 $\mathbb{F}_p$ tower sieve (density + wheels
  + nonexistence branch), 6.7 symbolic family intersection, 6.8 hardening protocol.
- **§8 drafted** as a density research plan: multidegree count ($E_4$ type $(2,2,2,4)$,
  exponent $H^{-2}$ yet families exist → special loci; $E_5$ type $(2^7,4^3,8)$,
  $\mathcal{O}(18)$, $H^{-18}$), the $r_2$ second-moment refinement (lift [TLT] $E_3$
  → $E_4$), a 5-step plan. **To confirm:** the exact $V_5$ multidegree/codim/canonical
  class (stated as research-plan item, not asserted).

**Advisor-level strategic notes (unchanged risks).**
- The **theorem inventory is still light**: Thm 4.1 (easy cascade of CMSV), Thm 5.3
  (one-paragraph), Prop 3.1 (reformulation of BMR sym-perfect). The one genuinely
  fiddly original piece — Prop 5.1 surjectivity — is now done. **Direction 5 (a new
  required divisor $C^E_{16}$ specific to the towered class) is the only clear shot at
  a genuinely new arithmetic theorem**; gated on reading CMSV §3 (their local families
  foreclose some improvements — but not necessarily the towered ones).
- **Literature pass (W0) is still blocking** the "first/new" framing: is the recursion
  in Gloden/Dorwart–Brown/Chernick? Are Bremner's two families Chernick/Gloden
  rediscovered? Identify before any "we introduce."
- **The campaign (§7) is a pilot, not a run.** Everything computational above is
  design, not results.

## 5d. Enumeration/implementation insights (→ `drafts/gpu_algorithms.md`)

Found while drafting the GPU algorithm spec; all verified numerically.

- **The whole $E_4$ search collapses to one invariant.** At the norm-$L_1$ level
  (not the root level), with $z=m+ni$, $N(z)=L_1$, define $\nu(z):=\mathrm{Re}(z^2)=m^2-n^2$.
  Then $\nu = r\cdot s$ for the twisted root pair $(r,s)=(m-n,m+n)$ (norm $2L_1$), and
  the $E_4$ condition is exactly $\nu_1^2+\nu_2^2=\nu_3^2+\nu_4^2$. $\nu^2$ is invariant
  under units **and** conjugation, so it is well defined on octant-canonized points.
  ⇒ **$E_4$ enumeration = a disjoint-pair-sum collision on $D$ integers.** Verified on
  the exemplar ($\nu_0 = 11810488155 = A\cdot B$; both sides $=1207464164080651264050$).
- **The pair key FACTORS (the real lever).** For a pair with $g=\gcd(z_1,z_2)$, $h=z_1/g$:
  $L_{2,d} = 2\,\mathrm{Re}(g^4)\mathrm{Re}(h^4)$, and $L_{2,d}=8L_1^2-8L_{2,b}$ (both
  verified; $g=165+467i$, $h=259-450i$). Since $L_1$ is fixed per norm, $L_{2,b}$ and
  $L_{2,d}$ are interchangeable keys — but only $L_{2,d}$ exposes that **the key is a
  product of two independent, unit- and conjugation-invariant factors**, one per side of
  a splitting of the prime multiset. This is the §5 "prism" in computational form.
  Whether it gives a sub-quadratic collision search is **open**, and is the highest-value
  algorithmic question for the champion regime.
- **Workload shape (exact sieve counts, not estimates).** Admissible $L_1\le B$
  (products of split primes, $D\ge4$): $B{=}10^6\to$ 11,083 norms / 99K pairs;
  $10^7\to$ 131,437 / 1.34M; $10^8\to$ 1,464,114 / 16.8M. Growth $\approx11\times$
  per decade (norms), $\approx12.5\times$ (pairs) ⇒ extrapolated $B{=}10^{12}$:
  $\sim2\times10^{10}$ norms, $\sim4\times10^{11}$ pair-sums. **The campaign is
  feasible**; the binding constraint is not FLOPs.
- **The $D$-distribution is top-heavy, which sets the kernel design.** At $10^8$:
  $D{=}4$ is 75.5% of norms, $D\in\{4,6,8\}$ is 96.4%. $D{=}4$ has exactly **3**
  disjoint pair-partitions ⇒ a register-only, hash-free, branch-free kernel covers
  three-quarters of all work. Disjoint-pair-pairs number $D(D{-}1)(D{-}2)(D{-}3)/8$,
  so direct comparison beats hashing up to $D\approx8$–$10$ (~97% of norms).
  ⇒ S-A is a **throughput problem over $10^{10}$ tiny uniform tasks**, not a
  compute-bound one.
- **64-bit is not enough.** The key $\nu_1^2+\nu_2^2 \le 2L_1^2$ is **71 bits already at
  the exemplar** and 81 bits at $B=10^{12}$. Policy: filter mod the Mersenne prime
  $2^{61}-1$, verify exactly. Crucially **modular filtering has no false negatives** —
  equal keys are always congruent — so the failure mode is extra work, never a silent
  miss. That is precisely the asymmetry an exclusion claim needs (§6.8).
- **Largest-prime-outer enumeration is the key structural choice.** Emit each norm once
  as $n = c\cdot p^e$ ($p$ = largest prime, $c$ from a small precomputed table of
  split-products $\le B/p$). This amortizes Cornacchia to **once per prime**
  ($3\times10^8$) instead of once per norm ($2\times10^{10}$), inherits $U(cp)$ from
  $U(c)$ with one Gaussian multiply per point, needs no factoring at all, and streams
  perfectly. The binding resource becomes the **factor-base table** ($\pi_1(1.5\times10^{10})
  \approx 3\times10^8$ primes $\approx$ 4 GB vs 8 GB device memory) — segment or compress it.
- **Hardware on hand:** RTX 4060 Laptop, 8 GB, sm_89, CUDA 12.9. `double` runs at 1/64
  rate ⇒ integer-only kernels. `code/native/` is a proven scaffold (GMP/CUDA/pybind/CTest)
  but a *sandbox*, not the enumerator.

### 5d.1 Second pass (rev. 2 of `gpu_algorithms.md`) — four new results

- **$\nu$ is INJECTIVE on canonized points** (proof: $m^2{+}n^2{=}L_1$ and
  $m^2{-}n^2{=}\nu$ determine $m,n$; checked on 300 random norms). Consequence:
  **two distinct pairs with equal key are automatically disjoint** — if they shared a
  point, $\nu_b^2=\nu_c^2 \Rightarrow b=c$. So *every* key collision is a genuine
  four-point configuration and the disjointness test disappears from every kernel.
  Also lets a point be stored as its $\nu$ alone (one word).
- **A4′ product form, verified generally at the norm-$L_1$ level** (12,088 pairs, 0
  mismatches): with $g$ = agree-part, $h$ = differ-part,
  $$\nu_1^2+\nu_2^2 = \mathrm{Re}(g^4)\,\mathrm{Re}(h^4) + L_1^2 .$$
  Derivation: $\mathrm{Re}(z^4)=2\nu^2-L_1^2$ and $z_1^4+z_2^4=g^4(h^4+\bar h^4)$.
  ⇒ **the $E_4$ condition is exactly a PRODUCT equality**
  $\mathrm{Re}(g_1^4)\mathrm{Re}(h_1^4)=\mathrm{Re}(g_2^4)\mathrm{Re}(h_2^4)$.
  *Status after exploration:* buys a ~2× constant (one multiply + two table loads vs.
  two squarings + add) and a **structural** prime-content filter (Thm 5.3 confines the
  primes to $\equiv\pm1\bmod8$, halving the basis; each $\mathrm{Re}(g^4)$ is a product
  of two $\mathbb{Z}[\sqrt2]$ norm values of size $\sim L_1^{1/2}$, cheap to factor).
  **No asymptotic gain yet** — $\sum_S 2^{s-2}=2^{2s-3}$ equals the naive pair count.
  Live direction: collisions are equalities of *factorizations over a restricted prime
  set*, so the natural framing is index-calculus (exponent vectors + linear algebra for
  4-term relations) rather than enumeration. Champion-regime tool only; in the
  exhaustive regime $D\le8$ and the naive key is already 3 ALU ops.
- **Modulus choice — plain $2^{64}$ wins; the Mersenne idea was wrong.** Key differences
  carry forced 2- and 3-adic content (Dilcher Rmk 3). Measured on 4,000 random
  differences: $\min v_2 = 6$, $\min v_3 = 2$, so $2^6\cdot3^2=576$ divides *every*
  difference. A mod-$2^k$ filter therefore loses exactly its bottom 6 bits:
  $2^{16}\to$ 10 effective bits (useless), $2^{32}\to$ 26 ($\sim6{\times}10^3$ false
  positives over the campaign, usable as a pre-filter), $2^{64}\to$ **58 effective bits,
  $\sim10^{-6}$ false positives, and free** (natural unsigned overflow). No-false-negative
  property is unaffected.
- **Multiplicities / imprimitivity (the subtle one).** For $L_1=p^bM$, a point with
  $p$-exponent $c$ is divisible by $p$ iff $1\le c\le b-1$, so the imprimitive share of
  $U$ is exactly $(b-1)/(b+1)$ — **1/3 at $b=2$** — and that family *is*
  $p\cdot U(M)$ (verified at $L_1=5^2\cdot13\cdot17$: $|U|=6$, imprimitive $=2=|U(13\cdot17)|$,
  sets coincide after scaling). The dedup rule is **not** "skip imprimitive points":
  - all four points imprimitive ⇒ $p^2\times$ a solution at $M$ ⇒ **duplicate, reject**
    (otherwise every solution reappears at $p^2M,q^2M,\dots$ and corrupts catalog counts
    and any density fit);
  - *some but not all* imprimitive ⇒ **primitive, keep** — this is precisely the S3 case
    (one pair uses $\mathfrak p\bar{\mathfrak p}=p$, another $\mathfrak p^2$), so
    discarding it reopens the hole Prop 5.1 exists to close.
  Implement as $\gcd(z_1,\dots,z_4)\ne1$ via a bitmask on the exponent vectors, on hits
  only. Non-squarefree norms are 17–21% of admissible norms but carry 35–37% of pair
  work; recoverable waste is only ~4%, so this is a **correctness fix, not an optimization**.
- **Scheduling by factor overlap is what makes it GPU-shaped** (Gil's suggestion). Order
  norms along the cofactor DFS tree, assign one block per subtree, keep the ancestor
  path's point sets in shared memory ($\le 8$ levels $\times D_{\max}{=}192 \times 8$B
  $= 12$ KB at $B{=}10^{12}$). A child's $U$ costs one Gaussian multiply per point from
  its parent, and the block emits thousands of norms with **zero global traffic**. This
  raises arithmetic intensity from 0.5 ops/byte (the naive "CPU generates, GPU meets"
  split, which is ~80× PCIe-bound and slower than not using the GPU) by ~3 orders of
  magnitude — it is what makes a device-resident design viable at all.
- $D_{\max}=192$ at $B=10^{12}$ (attained at $5^2\cdot13\cdot17\cdot29\cdot37\cdot41\cdot53\cdot61$).
  **$D$ stays small** — the single most important engineering fact for the kernels.

## 6. Verified computations (`drafts/verify_claims.py`)

- Exemplar $E_4$ (`A=252885,...,H=148453`): four pair-norms all $=66131993434$,
  $L_1=33065996717=13\cdot89\cdot173\cdot233\cdot709$. Recovered tower
  $L_2=489628056848329146064$, $L_3=175480010455650701584492675662518592000000$,
  $\sqrt{L_4}=40042900368028062136207226327668992000000$
  $=2^{14}3^65^67^211^2\cdot13\cdot17\cdot19\cdot29\cdot37^2\cdot43\cdot61\cdot107\cdot109\cdot139\cdot263\cdot389\cdot499$.
  **Full degree-16 identity $\prod(x^2-r_k^2)=(((x^2-L_1)^2-L_2)^2-L_3)^2-L_4$
  holds** (proof: both monic degree 16, identical).
- $L_{2,d}=-912832201535971887688$ (**negative** — proposal's sign typo confirmed)
  $=-2^3\cdot7\cdot17\cdot23^2\cdot31\cdot41\cdot47\cdot137\cdot401\cdot601\cdot919$;
  every odd prime $\equiv\pm1\pmod8$ ✓ (ℤ[√2] theorem).
- Ladder: $C'_8/2=151351200\mid\sqrt{L_4}$ ✓; $C'_4/2=18\mid\sqrt{L_3\pm\sqrt{L_4}}$,
  both perfect squares ✓.
- $E_5$: $\sqrt{L_5}\ge C'_{16}/2=2437428918743498865144960000$; review's
  $1.43\times10^{26}$ omitted the prime 17.
- **Ch.5 (§6 of the script):** exemplar has equal norms and $W = 0$ exactly
  (single-constraint, Prop 5.2); the surjectivity reconstruction ($\Phi^{-1}$ via the
  token argument) reproduces the inputs up to units for both the squarefree exemplar
  and a primitive $5^2\cdot13\cdot17$ quadruple with $q$-exponents $(2,1,0,2)$ — the
  $c=1$ entry being the pair that uses the rational factor $5 = \mathrm{q}\overline{\mathrm{q}}$
  (the non-squarefree/S3 case). The reconstruction code doubles as the §6 enumerator's
  non-squarefree recipe.

---

## 7. Reference index

**Obtained (in `sources/`):**
- `Bre08.pdf` — Bremner, Exp. Math. 17:4 (2008) 385–390. = ref [3]. The $E_4$ paper.
- `BMR09.pdf` — Borchert–McKenzie–Reinhardt, MFCS 2009. Gems, litter, sym-perfect.
- `CSMV23.pdf` — Coppersmith–Mossinghoff–Scheinerman–VanderKam, arXiv:2304.11254.
- `TLT26.pdf` — Tsai–Lee–Takahashi, arXiv:2606.07735. $E_3$ density.
- `Dil00.pdf` — Dilcher, Exp. Math. 9:3 (2000) 369–372. $E_3$/$E_4$ reductions;
  the refuted $E_4$-impossibility claim; $L_2$-difference $\mid 1152$.
- `RN19.pdf` — Raghavendran–Narayanan, *Mathematics* 7(3):227 (2019). PTE survey;
  history + bibliography spine for §1.
- `PTEProblemPaper2016Submitted.pdf` — **Filaseta & Markovich**, "Newton polygons
  and the PTE problem", *J. Number Theory* 174 (2017) 384–400. **2-adic values of
  the constant**: $\bar C_8 = 2^{e_1}3^35^27^2\cdot11\cdot13$ ($6\le e_1$),
  $\bar C_9 = 2^{e_2}3^{e_3}5^27^2\cdot11\cdot13\cdot17^{e_4}23^{e_5}29^{e_6}$ with
  $v_2(\bar C_9)=9$ proven; lists $\bar C_3..\bar C_7$. Writes **Tarry's 1913
  example exactly as our §2 counterexample**:
  $(x^2-5^2)(x^2-14^2)(x^2-23^2)(x^2-24^2)-(x^2-2^2)(x^2-16^2)(x^2-21^2)(x^2-25^2)
  =2^8\cdot3^3\cdot5^2\cdot7^2\cdot11\cdot13$ — corroborates §2 and pins the
  low-$n$ constants used in the §4 ladder. Relevant to §3/§5 2-power bookkeeping
  (their $k_1=k_1'$, $k_2\equiv k_2'\bmod4$ congruences).
- `borwein-sums1.pdf` — **Borwein & Choi**, "On Dirichlet series for sums of
  squares" (CECM preprint 2005; Ramanujan J.). Source of the **second-moment
  asymptotic $\sum_{n\le X}r_2(n)^2 = 4X\log X + 4\alpha_{BC}X + O(X^{2/3})$**
  that drives the §8 density (and [TLT]'s $E_3$ count).
- `PTE007.pdf` — **Wróblewski**, "A Collection of Numerical Solutions of Multigrade
  Equations Related to the PTE Problem" (v7, 2009), 68 pp. A **symmetric-multigrade
  solutions database**, indexed by max even exponent; focuses on sizes $\ge 10$
  (above our $E_4$ catalog at size 8, so overlap risk is low), confirms **no ideal
  solution beyond size 12** (matches [CMSV]). The §7 novelty-diligence table
  (alongside Chen Shuwen / Meyrignac). Its solution codes (e.g. `8.10.313`) are the
  Borwein–Lisoněk–Percival / Létac size-10 solutions. Cites Cipu, *Upper bounds for
  norms of products of binomials*, LMS J. Comput. Math. 7 (2004).

**To obtain (cited by the above, directly relevant; downloads were bot-walled):**
- **Borwein & Ingalls**, "The Prouhet–Tarry–Escott problem revisited",
  *L'Enseign. Math.* 40 (1994), 3–27. — the standard PTE survey; source of the
  degree-7 counterexample and the "ideal for every $n$?" conjecture. HIGH priority.
- **Dorwart & Brown**, "The Tarry–Escott problem", *Amer. Math. Monthly* 44 (1937),
  613–626. — equal power sums ⟺ polynomial difference degree (BMR Prop 3.4).
- **Chernick**, "Ideal solutions of the Tarry–Escott problem", *Amer. Math.
  Monthly* 44 (1937), 626–633. — classical symmetric constructions (check overlap
  with Bremner's families).
- **Rees & Smyth**, "On the constant in the Tarry–Escott problem", *Cinquante Ans
  de Polynômes*, LNM 1415 (1990), 196–208. — the required-divisor theory (CSMV
  Prop 2.1).
- **Caley**, Ph.D. thesis, "The Prouhet–Tarry–Escott problem", Waterloo 2012;
  and *Math. Comp.* 82 (2013) 1121–1137 (PTE over $\mathbb Z[i]$). — the ring case,
  overlaps our §5 Gaussian parameterization.
- **Choudhry & Wróblewski**, "Ideal solutions ... of degree eleven ...",
  *Hardy–Ramanujan J.* 31 (2008). — the size-12 families (CSMV eqn 7).
- **Existing solution tables** (further novelty diligence for the catalog, beyond
  Wróblewski `PTE007.pdf`): Chen Shuwen's multigrade site; Meyrignac's *eslp*
  database.
