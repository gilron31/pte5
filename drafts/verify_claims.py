import sympy as sp
from sympy import factorint, Integer, sqrt, Rational, Poly, symbols

y = symbols('y')

print("="*70)
print("1. EXEMPLAR E4 SOLUTION: tower recovery + polynomial identity")
print("="*70)
# 4 Gaussian integers (pairs) -> 8 positive roots
A,B,C,D,E,F,G,H = 252885,46703,195203,167415,249703,61485,209985,148453
coords = [A,B,C,D,E,F,G,H]
norms = [A*A+B*B, C*C+D*D, E*E+F*F, G*G+H*H]
print("pair norms:", norms, "all equal:", len(set(norms))==1)
twoL1 = norms[0]; L1 = twoL1//2
print("2L1 =", twoL1, " L1 =", L1, " L1 int:", twoL1%2==0)
print("L1 factor:", factorint(L1))

# Build Q(y) = prod (y - r_i^2), degree 8 in y = x^2
Q = sp.prod([ (y - r*r) for r in coords ])
Q = sp.Poly(sp.expand(Q), y)
print("Q monic deg", Q.degree())

def poly_sqrt_peel(P):
    """P monic deg 2d in y. Return (S, const) with S monic deg d, S^2 - P = const (should be constant)."""
    d = P.degree()//2
    coeffs = P.all_coeffs()  # high->low, length 2d+1
    # Build S = y^d + s_{d-1} y^{d-1}+...+s_0 by matching top d+1 coeffs of S^2 to P
    s = [Integer(1)] + [None]*d  # s[0]=lead
    # S^2 coeff of y^{2d-k} = sum_{i+j=k} s_i s_j ; solve for s_k
    for k in range(1, d+1):
        # coeff of y^{2d-k} in P
        target = coeffs[k]
        acc = Integer(0)
        for i in range(1, k):
            j = k-i
            if i<=d and j<=d:
                acc += s[i]*s[j]
        # term 2*s_0*s_k = 2*s_k
        sk = (target - acc)/2
        s[k]=sk
    S = sp.Poly([s[i] for i in range(d+1)], y)
    const = sp.expand((S*S - P).as_expr())
    return S, const

# Peel L4:  Q = S1^2 - L4  =>  L4 = S1^2 - Q  (the returned const)
S1, negL4 = poly_sqrt_peel(Q)
L4 = negL4
print("L4 constant residue is constant:", not negL4.free_symbols)
sqrtL4 = sp.sqrt(L4)
print("L4 =", L4)
print("sqrt(L4) =", sqrtL4, " perfect square:", sqrtL4.is_Integer)
# Peel L3: S1 = S2^2 - L3  =>  L3 = S2^2 - S1  (the returned const)
S2, negL3 = poly_sqrt_peel(S1)
L3 = negL3
print("L3 =", L3, " constant:", not negL3.free_symbols)
# S2 = (y-L1)^2 - L2 = y^2 -2L1 y + (L1^2 - L2)
s2c = S2.all_coeffs()
print("S2 coeffs:", s2c)
L1_chk = -s2c[1]/2
L2 = L1_chk**2 - s2c[2]
print("L1 from S2:", L1_chk, " matches:", L1_chk==L1)
print("L2 =", L2)

print("\n-- compare to review's recovered tower --")
print("L2 review: 489628056848329146064  match:", L2==489628056848329146064)
print("L3 review: 175480010455650701584492675662518592000000  match:", L3==175480010455650701584492675662518592000000)
print("sqrtL4 review: 40042900368028062136207226327668992000000  match:", sqrtL4==40042900368028062136207226327668992000000)

# Full polynomial identity check: prod(x^2-r^2) == (((x^2-L1)^2-L2)^2-L3)^2-L4
x = symbols('x')
lhs = sp.expand(sp.prod([ (x*x - r*r) for r in coords ]))
rhs = sp.expand((((x**2-L1)**2-L2)**2-L3)**2-L4)
print("FULL POLY IDENTITY holds:", sp.expand(lhs-rhs)==0)

print("\n-- sqrt(L4) factorization --")
print(factorint(sqrtL4))

print("\n"+"="*70)
print("2. L_{2,d} sign + Z[sqrt2] parity")
print("="*70)
# L2d(a+ib,c+id) = Re((a+ib)^4+(c+id)^4). Use pair (A,B) and (C,D) for group1
def L2d(a,b,c,d):
    return sp.re((a+b*sp.I)**4) + sp.re((c+d*sp.I)**4)
val = L2d(A,B,C,D)
print("L2d(group1) =", val)
val2 = L2d(E,F,G,H)
print("L2d(group2) =", val2, " equal-magnitude:", abs(val)==abs(val2))
print("sign negative:", val<0)
fac = factorint(abs(val))
print("|L2d| factor:", fac)
print("odd primes = +-1 mod 8:", all((p%8 in (1,7)) for p in fac if p!=2))

print("\n"+"="*70)
print("3. DIVISIBILITY LADDER on exemplar (E4)")
print("="*70)
Cp8_over2 = 2**5 * 3**3 * 5**2 * 7**2 * 11 * 13
print("C'_8/2 =", Cp8_over2)
print("C'_8/2 | sqrt(L4):", sqrtL4 % Cp8_over2 == 0)
# level-3: C'_4/2 = 18 | sqrt(L3 +- sqrt(L4))
Cp4_over2 = (2**2*3**2)//2
print("C'_4/2 =", Cp4_over2)
for sgn in (+1,-1):
    inner = L3 + sgn*sqrtL4
    r = sp.sqrt(inner)
    print(f"  sqrt(L3 {'+' if sgn>0 else '-'} sqrtL4) integer:", r.is_Integer, " /18 int:", (r % 18==0) if r.is_Integer else 'NA')

print("\n"+"="*70)
print("4. E5 LOWER BOUNDS from C'_16 (Table 1, CSMV)")
print("="*70)
# Table 1 n=16: C_16 divisors, plus additional for C'_16
C16_full   = 2**11 * 3**6 * 5**4 * 7**3 * 11**2 * 13**2 * 17 * 19 * 23   # C_16 (incl boxed 23)
add_Cp16   = 29 * 37 * 41 * 43 * 53
Cp16 = C16_full * add_Cp16
Cp16_no17 = (C16_full//17) * add_Cp16   # review's version (dropped 17)
print("C'_16          =", Cp16, " ~ %.3e"%float(Cp16))
print("C'_16/2        =", Cp16//2, " ~ %.3e"%float(Cp16//2))
print("C'_16(no17)/2  ~ %.3e"%(float(Cp16_no17)/2), " (review claimed 1.43e26)")
sqrtL5_min = Rational(Cp16,2)
print("\nHEADLINE: sqrt(L5) >= C'_16/2 ~ %.3e"%float(sqrtL5_min))
print("         => L5 >= (C'_16/2)^2 ~ %.3e"%float(sqrtL5_min**2))
# derived L1 from ladder magnitude chain sqrt(L5) <= L1^8  => L1 >= (sqrtL5)^(1/8)
L1_ladder = float(sqrtL5_min)**(1/8)
print("ladder-only L1 >= (C'16/2)^(1/8) ~ %.3e (weak)"%L1_ladder)
# from CMSV exhaustive size-16 search: no sym solution height<=850 => any E5 has r_max>850
print("CMSV height: any E5 -> size16 sym sol, height r_max>850 => L1 >= 850^2/2 =", 850**2//2)
print("  (r_max=sqrt(2L1) approx). This BEATS ladder-only L1 bound.")

print("\n"+"="*70)
print("5. TOWER RECURSION (Prop 3.1) on exemplar E4 -> E3")
print("="*70)
xx = symbols('x')
uu = [(A*A-B*B)//2, (C*C-D*D)//2, (E*E-F*F)//2, (G*G-H*H)//2]
print("u =", uu)
print("u0^2+u1^2 == u2^2+u3^2 == 2L2 :",
      uu[0]**2+uu[1]**2 == 2*L2 and uu[2]**2+uu[3]**2 == 2*L2)
lhs_r = sp.expand(sp.prod([(xx**2 - uk**2) for uk in uu]))
rhs_r = sp.expand(((xx**2 - L2)**2 - L3)**2 - L4)   # L4 here is the true L4 (perfect square)
print("prod(x^2-u_k^2) == ((x^2-L2)^2-L3)^2 - L4 (derived E3):", sp.expand(lhs_r-rhs_r)==0)
print("congruum L1 +- u_k all perfect squares:",
      all(sp.sqrt(L1 + s*uk).is_Integer for uk in uu for s in (1,-1)))

print("\n"+"="*70)
print("6. CH.5: single-constraint W=0 + SURJECTIVITY of the parameterization")
print("="*70)
import math

# --- minimal exact Gaussian-integer arithmetic (re,im) tuples ---
def gmul(z, w): return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])
def gconj(z):   return (z[0], -z[1])
def gnorm(z):   return z[0]*z[0]+z[1]*z[1]
def gdiv(z, w):
    """z/w if exact in Z[i], else None."""
    n = gnorm(w)
    a = z[0]*w[0]+z[1]*w[1]      # Re(z*conj(w))
    b = -z[0]*w[1]+z[1]*w[0]     # Im(z*conj(w))
    if a % n == 0 and b % n == 0: return (a//n, b//n)
    return None
def gval(z, q):
    """(valuation of q in z, cofactor)."""
    v = 0
    while True:
        d = gdiv(z, q)
        if d is None: return v, z
        z = d; v += 1
def is_unit(z):      return gnorm(z) == 1
def associate(z, w): 
    d = gdiv(z, w)
    return d is not None and is_unit(d)
def split_prime(p):
    """q=(a,b), a>=b>0, a^2+b^2=p, for p==1 mod 4."""
    for a in range(1, math.isqrt(p)+1):
        b2 = p - a*a
        b = math.isqrt(b2)
        if b > 0 and b*b == b2:
            return (max(a, b), min(a, b))
    raise ValueError(p)

# the fixed conjugation pattern of the paper's quadruplet parameterization:
#   n1 conjugates none; n2 -> X4X5X6X7; n3 -> X2X3X6X7; n4 -> X1X3X5X7
CONJ = {2: {4,5,6,7}, 3: {2,3,6,7}, 4: {1,3,5,7}}
def phi(bk):
    """Phi(X_0..X_7) -> (n1,n2,n3,n4)."""
    def build(cset):
        z = (1,0)
        for t in range(8):
            xt = bk[t]
            z = gmul(z, gconj(xt) if t in cset else xt)
        return z
    return (build(set()), build(CONJ[2]), build(CONJ[3]), build(CONJ[4]))

def bucket_of(sigma):
    """subset sigma of {2,3,4} -> bucket index t (n2<->4, n3<->2, n4<->1)."""
    return 4*(2 in sigma) + 2*(3 in sigma) + 1*(4 in sigma)

def reconstruct(quad):
    """Given 4 equal-norm Gaussian ints, build X_0..X_7 with Phi ~ quad (Prop 5.1)."""
    n1 = quad[0]
    N = gnorm(n1)
    assert all(gnorm(z) == N for z in quad), "not equal-norm"
    bk = [(1,0)]*8
    for p, e in factorint(N).items():
        if p == 2:                                  # ramified: (1+i)^a common -> X_0
            a = gval(n1, (1,1))[0]
            for _ in range(a): bk[0] = gmul(bk[0], (1,1))
        elif p % 4 == 3:                            # inert: p^(e/2) common -> X_0
            for _ in range(e//2): bk[0] = gmul(bk[0], (p,0))
        else:                                       # split: token-distribute q,qbar
            q = split_prime(p); b = e
            c = [gval(z, q)[0] for z in quad]        # (c1,c2,c3,c4) in {0..b}^4
            for r in range(b):                       # r-th token
                m = [1 if r < c[j] else 0 for j in range(4)]   # m[0]=c1-bit ...
                if m[0] == 1:                        # q-token
                    sigma = {j for j in (2,3,4) if m[j-1] == 0}
                    fac = q
                else:                                # qbar-token
                    sigma = {j for j in (2,3,4) if m[j-1] == 1}
                    fac = gconj(q)
                t = bucket_of(sigma)
                bk[t] = gmul(bk[t], fac)
    return bk

def W(quad):
    (A,B),(C,D),(E,Fx),(G,Hh) = quad
    return A*A*B*B + C*C*D*D - E*E*Fx*Fx - G*G*Hh*Hh

# 6a. single constraint on the exemplar (equal norms + W=0)
exemplar = [(A,B),(C,D),(E,F),(G,H)]
print("exemplar equal norms:", len({gnorm(z) for z in exemplar})==1,
      " W = A^2B^2+C^2D^2-E^2F^2-G^2H^2 =", W(exemplar), " (==0:", W(exemplar)==0, ")")

# 6b. surjectivity: reconstruct buckets, check Phi ~ inputs up to units
def check_surj(quad, label):
    bk = reconstruct(quad)
    out = phi(bk)
    ok = all(associate(quad[j], out[j]) for j in range(4))
    occ = [t for t in range(8) if not is_unit(bk[t])]
    print(f"  [{label}] Phi(X) ~ inputs (up to units): {ok}   occupied buckets: {occ}")
    return ok

print("SURJECTIVITY reconstruction (Prop 5.1):")
ok1 = check_surj(exemplar, "exemplar (squarefree norm 2*13*89*173*233*709)")

# non-squarefree test: norm 5^2 * 13 * 17 = 5525, with a root pair using the REAL
# factor 5 = q*qbar (c=1) and another using q^2 (c=2) -> the primitive p^2 case (S3).
nsq = [(-14,73), (70,25), (22,-71), (74,7)]   # built from q=2+i,(3+2i),(4+i); see notes
print("  non-squarefree norms all == 5525:", {gnorm(z) for z in nsq})
q5 = split_prime(5)
print("  q=2+i exponents (c1..c4) incl the c=1 'real 5' pair:",
      [gval(z, q5)[0] for z in nsq])
ok2 = check_surj(nsq, "non-squarefree 5^2*13*17 (S3 primitive case)")

print("\nALL CH.5 CHECKS PASS:", (W(exemplar)==0) and ok1 and ok2)
