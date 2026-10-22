"""G1.2: the transversal type of V_5's singularities along each stratum of N12.

Near a point over (1 : d) whose zero leaves are r_{2k+1} = 0 (k in K), V_5 is
locally (etale in all other square roots) the complete intersection
    X_K = { V_4(d) = 0,  y_k^2 = (1 - d_k)/2  (k in K) }  in  C^{8+|K|}.
We compute its singular scheme (maximal minors of the Jacobian), keep only the
components through y_K = 0 (the stratum) and away from d_k = +-1 (k not in K),
and report for each component its dimension and whether the singular scheme is
REDUCED there. Reduced along a curve means the transversal 3-fold singularity
has Tjurina number 1, i.e. is an ordinary double point, which is canonical.
Usage: python3 v5_local_type.py
"""
import subprocess

V4 = ['d0^2+d1^2-d2^2-d3^2', 'd0^2+d1^2-d4^2-d5^2', 'd0^2+d1^2-d6^2-d7^2',
      'd0^2*d1^2+d2^2*d3^2-d4^2*d5^2-d6^2*d7^2']
STRATA = [(0,), (0, 2), (0, 4), (0, 1, 2, 3), (0, 1, 4, 5), (0, 2, 4, 6)]


def script(K):
    ys = [f'y{k}' for k in K]
    others = [k for k in range(8) if k not in K]
    eqs = V4 + [f'2*y{k}^2 - (1 - d{k})' for k in K]
    n = len(eqs)
    return f'''LIB "primdec.lib"; LIB "elim.lib";
option(redSB);
ring r = 0, (d0,d1,d2,d3,d4,d5,d6,d7,{",".join(ys)}), dp;
ideal I = {", ".join(eqs)};
ideal S = I + minor(jacob(I), {n});
// on the stratum: y_K = 0 locally; away from the other branch divisors and from d = 0
poly opn = {"*".join(f"(1 - d{k}^2)" for k in others) or "1"};
def S0 = sat(S, ideal(opn));
ideal T; if (typeof(S0) == "list") {{ T = S0[1]; }} else {{ T = S0; }}
T = std(T + ideal({", ".join(ys)}) * 0);
list P = primdecGTZ(T);
int i;
for (i = 1; i <= size(P); i++) {{
  ideal Q = std(P[i][1]); ideal R = std(P[i][2]);
  int onstratum = (reduce({ys[0]}, R) == 0);
  if (onstratum) {{
    int red = (size(reduce(R, Q)) == 0);
    print("COMPONENT dim " + string(dim(R)) + " degree " + string(mult(R)) + " reduced " + string(red));
  }}
  kill onstratum;
}}
quit;
'''


if __name__ == '__main__':
    for K in STRATA:
        open('v5_local_tmp.sing', 'w').write(script(K))
        try:
            out = subprocess.run(['Singular', '-q', 'v5_local_tmp.sing'], capture_output=True, text=True, timeout=1800).stdout
        except subprocess.TimeoutExpired:
            print(f'K = {K}: timeout'); continue
        comps = [l.strip() for l in out.splitlines() if l.startswith('COMPONENT')]
        print(f'K = {K!s:<16} components through the stratum: {len(comps)}')
        for c in comps:
            print('    ' + c)
