\\ Scan phi-fibers on Bremner's slice for rational points with a different E_4.
\\ phi(s:t:u) = (L1:L2:L3); for P = (1:t0:u0) find all rational (1:t:u) with
\\ phi equal, via the resultant in u, and report those whose j4 = S4^2/L1^8
\\ differs from P's (a different E_4 over Qbar with the same L1, L2, L3).
\\ Run: gp -q phi_scan.gp
inv(S, T, U) =
{
  my(A = -S + T, C = S + T, E = -S + U, G = S + U);
  my(b = 5*S^2 + 2*S*T + U^2, d = 5*S^2 - 2*S*T + U^2);
  my(f = 5*S^2 + 2*S*U + T^2, h = 5*S^2 - 2*S*U + T^2);
  my(L1 = A^2 + b, L2 = (A^2 - b)^2 + (C^2 - d)^2);
  my(W1 = A^2*b - C^2*d, W2 = E^2*f - G^2*h);
  [L1, L2, W1^2 + W2^2, W1^2 - W2^2];
}
INV = inv(1, 't, 'u);
partners(t0, u0) =
{
  my(a = inv(1, t0, u0), out = List());
  if (a[1] == 0, return(out));
  my(e1 = INV[2]*a[1]^2 - a[2]*INV[1]^2, e2 = INV[3]*a[1]^4 - a[3]*INV[1]^4);
  my(R = polresultant(e1, e2, 'u), j4P = a[4]^2 / a[1]^8, F);
  if (R == 0, return(out));
  F = factor(R)[,1];
  for (i = 1, #F,
    if (poldegree(F[i], 't) != 1, next);
    my(tt = -polcoef(F[i], 0, 't) / polcoef(F[i], 1, 't));
    my(g = gcd(subst(e1, 't, tt), subst(e2, 't, tt)));
    if (poldegree(g, 'u) < 1, next);
    my(G = factor(g)[,1]);
    for (k = 1, #G,
      if (poldegree(G[k], 'u) != 1, next);
      my(uu = -polcoef(G[k], 0, 'u) / polcoef(G[k], 1, 'u), q = inv(1, tt, uu));
      if (q[1] == 0, next);
      my(j4 = q[4]^2 / q[1]^8);
      if (j4 != j4P, listput(out, [tt, uu, j4 / j4P]))));
  out;
}
{
my(N = 25, nhit = 0, nspecial = 0, t0 = getwalltime());
for (t0_ = -N, N, for (u0 = 1, N,
  if (gcd(t0_, u0) != 1 && t0_ != 0, next);
  my(p = partners(t0_, u0));
  if (#p, nhit++; printf("P = (1:%d:%d): %d rational partners with a different E_4, e.g. %s\n", t0_, u0, #p, p[1]))));
printf("scanned |t| <= %d, 1 <= u <= %d: %d base points with a nontrivial rational partner (%d ms)\n",
       N, N, nhit, getwalltime() - t0);
}
