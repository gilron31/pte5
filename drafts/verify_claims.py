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
