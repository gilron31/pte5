default(parisizemax, 2^31);
\\ Certify rank J(C_S)(Q) = 1 for genus-2 quotients C_S: y^2 = prod_{e in S} (t + e)
\\ of the bottom-lift curve over Symes' E_4:
\\   (1) the descent bound is 1 with every local image at full size;
\\   (2) J(Q)_tors = J[2]: gcd of #J(F_p) over good primes, odd part trivial
\\       and 2-part 16 (all 2-torsion is rational, so 16 | #J(Q)_tors);
\\   (3) the trivial-lift point P0 = (t = oo) gives a class whose x - T image is not
\\       in the span of the 2-torsion images, so it is not torsion.
\\ Then rank = 1 exactly.  Output: one line per certified subset.
read("g2descent.gp");
R = [131, 343, 11, 367, 101, 353, 77, 359];
signed = concat(R, -R);

certify(S) =
{
  my(e = vector(6, i, -signed[S[i]]), g = g2rank(e), a = odd_model(e), f, h, N = 0);
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
{
my(cands = [[1, 5, 6, 11, 12, 14], [6, 7, 8, 9, 14, 15], [1, 9, 11, 12, 14, 15], [2, 3, 8, 9, 11, 13],
             [3, 5, 8, 12, 13, 16], [2, 3, 4, 7, 11, 14], [5, 8, 9, 11, 14, 15], [3, 6, 7, 8, 9, 13]], ok = 0);
foreach(cands, S,
  my(r = certify(S));
  printf("S = %s (roots %s): %s\n", S, vector(6, i, -signed[S[i]]), if (r[1], "rank J(Q) = 1 CERTIFIED", r[2]));
  ok += r[1]);
printf("%d of %d candidates certified rank 1\n", ok, #cands);
}
