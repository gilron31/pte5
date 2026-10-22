fn = "fam1_F_gp.txt";
\\ Rational singular points of a family's image curve F(J1, J2) = 0, via resultants.
\\ Usage: gp -q sing_rational.gp  (reads FAMILY_F_gp.txt via the variable fn)
default(parisizemax, 2^33);
F = eval(Str(read(fn)));
{
my(t0 = getwalltime(), FA = deriv(F, 'J1), FB = deriv(F, 'J2), R1, R2, g, fa, pts = List());
R1 = polresultant(F, FB, 'J2);
R2 = polresultant(FA, FB, 'J2);
g = gcd(R1, R2);
printf("deg R1 %d, deg R2 %d, deg gcd %d (%d s)\n", poldegree(R1), poldegree(R2), poldegree(g), (getwalltime() - t0) \ 1000);
fa = factor(g)[,1];
for (i = 1, #fa,
  if (poldegree(fa[i], 'J1) != 1, next);
  my(j1 = -polcoef(fa[i], 0, 'J1) / polcoef(fa[i], 1, 'J1));
  my(h = gcd(gcd(subst(F, 'J1, j1), subst(FA, 'J1, j1)), subst(FB, 'J1, j1)));
  if (poldegree(h, 'J2) < 1, next);
  my(H = factor(h)[,1]);
  for (k = 1, #H, if (poldegree(H[k], 'J2) == 1, listput(pts, [j1, -polcoef(H[k], 0, 'J2) / polcoef(H[k], 1, 'J2)]))));
printf("rational singular points: %d\n", #pts);
for (i = 1, #pts, print("  ", pts[i]));
}
