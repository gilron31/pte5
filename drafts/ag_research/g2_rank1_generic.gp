\\ Generic rank-1 certificate for y^2 = prod_{i=1..6} (t - e_i), e_i distinct rationals:
\\ returns [1, dim Sel] if rank J(Q) = 1 is proven (descent bound 1 with full local
\\ images, torsion = J[2], trivial-lift point of infinite order), else [0, reason].
\\ See g2_rank1_certify.gp for the argument.
read("g2descent.gp");
certify_e(e) =
{
  my(g = g2rank(e), a = odd_model(e), f, h, N = 0);
  if (g[1] != 1, return([0, Str("bound ", g[1])]));
  foreach(g[4], ld, if (ld[2] < ld[3], return([0, Str("local image short at ", ld[1])])));
  \\ (2) torsion: gcd of #J(F_p) over good primes
  f = prod(i = 1, 6, 't - e[i]);
  forprime(p = 3, 400,
    if (vecsearch(g[3], p), next);
    my(fp = Mod(f, p));
    if (poldisc(lift(fp)) % p == 0, next);
    my(cp = hyperellcharpoly(fp));
    N = gcd(N, subst(cp, 'x, 1)));
  my(oddpart = N >> valuation(N, 2));
  if (oddpart != 1 || valuation(N, 2) != 4, return([0, Str("torsion bound ", N)]));
  \\ (3) the point at t = oo is X = 0 in the odd model (W^2 = prod(-a_i) must be a square)
  my(y2 = prod(i = 1, 5, -a[i]));
  if (!issquare(y2), return([0, "X = 0 not a rational point?"]));
  \\ global image of P0 - oo: square classes of (0 - a_i) over <-1, S>
  my(gl = concat([-1], g[3]), vec(x) = concat([if (x < 0, 1, 0)], vector(#gl - 1, k, valuation(x, gl[k + 1]) % 2)));
  my(T = Mod(Mat(vector(5, j, concat(vector(5, i, vec(mu_pt(a, a[j])[i])))~)), 2));
  my(P = Mod(concat(vector(5, i, vec(-a[i])))~, 2));
  if (matrank(matconcat([T, P])) == matrank(T), return([0, "P0 image in 2-torsion span"]));
  [1, g[2]];
}
