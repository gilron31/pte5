\\ The E2 quotients whose rank ellrank left at [0,2] at default effort: order of the
\\ trivial-lift point P_inf = (0, c), root number, ellrank at effort 3; then rank tallies.
\\ Run: gp -q lift_triples_check.gp
default(parisizemax, 2^31);
read("lift_lib.gp");
R = [131, 343, 11, 367, 101, 353, 77, 359];
{
foreach([[1,4,6],[1,6,8],[2,3,4],[4,6,8],[4,5,6]], v,
  my(a = R[v[1]], b = R[v[2]], c = R[v[3]], P, e2, cc, Pinf, rk);
  P = prod(i = 2, 4, 1 + (0 - [0, a^2, b^2, c^2][i]) * 'x);
  cc = polcoef(P, 3);
  e2 = mkE2(a, b, c);
  Pinf = [0, cc];
  rk = ellrank(e2, 3);
  printf("%s  on curve %d  ord(P_inf) %d  root number %d  ellrank(effort 3) [%d,%d]\n",
         v, ellisoncurve(e2, Pinf), ellorder(e2, Pinf), ellrootno(e2), rk[1], rk[2]));
}
{
\\ tallies over all 56 triples
my(t1 = vector(6), t2 = vector(6), und = 0);
forsubset([8, 3], v,
  my(a = R[v[1]], b = R[v[2]], c = R[v[3]], r1 = ellrank(mkE1(a, b, c)), r2 = ellrank(mkE2(a, b, c)));
  if (r1[1] == r1[2], t1[r1[1] + 1]++);
  if (r2[1] == r2[2], t2[r2[1] + 1]++, und++));
print("E1 rank tallies (rank 0..5): ", t1);
print("E2 rank tallies (rank 0..5): ", t2, "  undetermined: ", und);
}
{
\\ the trivial lift mu = 0 on every E2: P_inf = (0, c) has infinite order (ellorder = 0)
my(fin = 0);
forsubset([8, 3], v,
  my(a = R[v[1]], b = R[v[2]], c = R[v[3]], P = prod(i = 2, 4, 1 - [0, a^2, b^2, c^2][i] * 'x), e2 = mkE2(a, b, c));
  if (ellorder(e2, [0, polcoef(P, 3)]) != 0, fin++));
print("E2 curves on which the trivial-lift point is torsion: ", fin, " of 56");
}
