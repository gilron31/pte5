default(parisizemax, 2^31);
read("g2_rank1_generic.gp");
\\ A3 calibration seed: Symes' derived E_3 (roots below), known lift t0 = 67405/12.
rho = [4187, 5607, 4767, 5123];
signed = concat(rho, -rho);
t0 = 67405/12;
{
forsubset([8, 6], v,
  my(e = vector(6, i, -signed[v[i]]));
  if (#setintersect(Set(e), Set(-e)) == 6, next);
  my(g = g2rank(e), y2 = prod(i = 1, 6, t0 - e[i]));
  printf("%s: rank bound %d (local dims ok: %d); t0 on curve: %d\n", Vec(v), g[1],
         #select(ld -> ld[2] < ld[3], g[4]) == 0, issquare(y2)));
}
