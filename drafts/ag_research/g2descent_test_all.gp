default(parisizemax, 2^31);
\\ Validation of g2descent.gp on all 56 split Jacobians (symmetric S = {+-a,+-b,+-c}
\\ over Symes' roots): the bound must be >= rank E1 + rank E2 (when those are proven).
read("g2descent.gp");
read("../ag_pilots/lift_lib.gp");
R = [131, 343, 11, 367, 101, 353, 77, 359];
{
my(n = 0, sharp = 0, viol = 0, undec = 0, short = 0);
forsubset([8, 3], v,
  my(a = R[v[1]], b = R[v[2]], c = R[v[3]], e = [-a, a, -b, b, -c, c]);
  my(r1 = ellrank(mkE1(a, b, c), 3), r2 = ellrank(mkE2(a, b, c), 3), g = g2rank(e));
  if (#select(ld -> ld[2] < ld[3], g[4]), short++; next);
  n++;
  if (g[1] < r1[1] + r2[1], viol++; printf("VIOLATION %s: bound %d < proven lower bound %d\n", v, g[1], r1[1] + r2[1]));
  if (r1[1] == r1[2] && r2[1] == r2[2],
    if (g[1] == r1[1] + r2[1], sharp++),
    undec++));
printf("56 split Jacobians: %d usable (%d with short local image skipped); violations %d; descent sharp in %d of %d with known rank; %d with undetermined elliptic ranks\n",
       n, short, viol, sharp, n - undec, undec);
}
