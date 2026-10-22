\\ Sanity check for phi_scan.gp: all rational fiber points, no j4 filter.
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
allpts(t0, u0) =
{
  my(a = inv(1, t0, u0), n = 0, e1, e2, R, F, tt, g, G, uu);
  e1 = INV[2]*a[1]^2 - a[2]*INV[1]^2;
  e2 = INV[3]*a[1]^4 - a[3]*INV[1]^4;
  R = polresultant(e1, e2, 'u);
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
      if (inv(1, tt, uu)[1] != 0, n++)));
  [poldegree(R), n];
}
print("(1:3:7): [deg resultant, rational fiber points] = ", allpts(3, 7));
print("(1:2:5): ", allpts(2, 5));
print("(1:5:13): ", allpts(5, 13));
