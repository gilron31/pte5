\\ Independent check of N8c: intersect the two family image curves in the (J1,J2)-plane.
\\ Run: gp -q cross_images.gp
default(parisizemax, 2^32);
F1 = eval(Str(read("fam1_F_gp.txt")));
F2 = eval(Str(read("fam2_F_gp.txt")));
{
my(t0 = getwalltime(), R = polresultant(F1, F2, 'J2), fa, cand = List());
printf("resultant in J1: degree %d (%d s)\n", poldegree(R, 'J1), (getwalltime() - t0) \ 1000);
fa = factor(R)[,1];
printf("factor degrees: %s\n", vecsort(vector(#fa, i, poldegree(fa[i], 'J1))));
for (i = 1, #fa,
  if (poldegree(fa[i], 'J1) != 1, next);
  my(j1 = -polcoef(fa[i], 0, 'J1) / polcoef(fa[i], 1, 'J1));
  my(g = gcd(subst(F1, 'J1, j1), subst(F2, 'J1, j1)));
  if (poldegree(g, 'J2) < 1, next);
  my(G = factor(g)[,1]);
  for (k = 1, #G, if (poldegree(G[k], 'J2) == 1, listput(cand, [j1, -polcoef(G[k], 0, 'J2) / polcoef(G[k], 1, 'J2)]))));
printf("rational intersection points of the two image curves: %d\n", #cand);
for (i = 1, #cand, print("  ", cand[i]));
}
