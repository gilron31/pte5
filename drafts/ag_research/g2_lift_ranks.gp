default(parisizemax, 2^31);
\\ G5: proven rank bounds for genus-2 quotients y^2 = prod_{e in S}(t + e) of the
\\ bottom-lift curve over Symes' E_4, S a 6-subset of the 16 signed roots.
\\ A bound of 1 means rank exactly 1 (the trivial lift mu = 0 gives a point of
\\ infinite order -- checked separately), i.e. Chabauty-ready.
read("g2descent.gp");
R = [131, 343, 11, 367, 101, 353, 77, 359];
signed = concat(R, -R);
{
my(tally = vector(10), ones = List(), n = 0, t0 = getwalltime());
setrand(1);
while (n < 400,
  my(v = vecsort(vector(6, i, 0)), idx = List());
  \\ random 6-subset of 1..16, rejecting symmetric ones (closed under e -> -e)
  my(perm = numtoperm(16, random(16!)), S = vecsort(vector(6, i, perm[i])));
  my(e = vector(6, i, -signed[S[i]]));
  my(sym = (#setintersect(Set(e), Set(-e)) == 6));
  if (sym, next);
  n++;
  my(g = g2rank(e), b = g[1]);
  if (b >= 0 && b < 10, tally[b + 1]++);
  if (b == 1, listput(ones, S)));
printf("400 random non-symmetric 6-subsets: rank-bound tally (bound 0..9): %s  (%d s)\n", tally, (getwalltime() - t0) \ 1000);
printf("first subsets with proven rank <= 1: %s\n", vector(min(8, #ones), i, ones[i]));
}
