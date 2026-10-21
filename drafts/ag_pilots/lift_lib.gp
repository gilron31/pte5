\\ Bottom-lift curves for an E_4 seed with positive roots R.
\\ quartic_E(al): Jacobian (Weierstrass) of y^2 = prod (t - al[i]), al[1] sent to infinity.
quartic_E(al) =
{
  my(a1 = al[1], c, P, s = 'x);
  P = prod(i = 2, 4, 1 + (a1 - al[i]) * s);
  c = polcoef(P, 3);
  ellinit([0, polcoef(P, 2), 0, polcoef(P, 1) * c, c^2]);
}
mkE1(a, b, c) =
{
  my(P = (x - a^2) * (x - b^2) * (x - c^2));
  ellinit([0, polcoef(P, 2), 0, polcoef(P, 1), polcoef(P, 0)]);
}
mkE2(a, b, c) = quartic_E([0, a^2, b^2, c^2]);
