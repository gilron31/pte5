\\ Wider scan of phi-fibers, plus the known slice points themselves.
read("phi_scan_check.gp");
partners(t0, u0) =
{
  my(a = inv(1, t0, u0), out = List(), nrat = 0, e1, e2, R, F, tt, g, G, uu, q, j4P);
  if (a[1] == 0, return([out, 0]));
  j4P = a[4]^2 / a[1]^8;
  e1 = INV[2]*a[1]^2 - a[2]*INV[1]^2;
  e2 = INV[3]*a[1]^4 - a[3]*INV[1]^4;
  R = polresultant(e1, e2, 'u);
  if (R == 0, return([out, -1]));
  F = factor(R)[,1];
  for (i = 1, #F,
    if (poldegree(F[i], 't) != 1, next);
    tt = -polcoef(F[i], 0, 't) / polcoef(F[i], 1, 't);
    g = gcd(subst(e1, 't, tt), subst(e2, 't, tt));
    if (poldegree(g, 'u) < 1, next);
    G = factor(g)[,1];
    for (k = 1, #G,
      if (poldegree(G[k], 'u) != 1, next);
      uu = -polcoef(G[k], 0, 'u) / polcoef(G[k], 1, 'u);
      q = inv(1, tt, uu);
      if (q[1] == 0, next);
      nrat++;
      if (q[4]^2 / q[1]^8 != j4P, listput(out, [tt, uu]))));
  [out, nrat];
}
{
my(N = 400, hits = 0, special = List(), t0 = getwalltime());
for (tt0 = -N, N, for (uu0 = 1, N,
  if (gcd(tt0, uu0) != 1, next);
  my(r = partners(tt0, uu0));
  if (#r[1], hits++; printf("NONTRIVIAL: (1:%d:%d) partners %s\n", tt0, uu0, r[1]));
  if (r[2] > 8, listput(special, [tt0, uu0, r[2]]))));
printf("box |t| <= %d, 1 <= u <= %d: %d base points with a different-E_4 rational partner; %d special (more than 8 rational, same E_4); %d s\n",
       N, N, hits, #special, (getwalltime() - t0) \ 1000);
print("first special points: ", vector(min(12, #special), i, special[i]));
}
\\ the slice points (s:t:u) of the known slice E_4 (from stu_images.py), in the chart s = 1
{
my(K = [[74250, 135735, 120953]]);
foreach(K, v, my(r = partners(v[2] / v[1], v[3] / v[1]));
  printf("exemplar image (%d:%d:%d): %d rational fiber points, nontrivial partners %s\n", v[1], v[2], v[3], r[2], r[1]));
}
