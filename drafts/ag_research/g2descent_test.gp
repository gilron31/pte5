default(parisizemax, 2^31);
\\ Test the genus-2 descent on symmetric subsets S = {+-a, +-b, +-c}, where
\\ J ~ E1 x E2 with E1: Y^2 = (X-a^2)(X-b^2)(X-c^2), E2: v^2 = u(u-a^2)(u-b^2)(u-c^2).
read("g2descent.gp");
read("../ag_pilots/lift_lib.gp");
R = [131, 343, 11, 367, 101, 353, 77, 359];
{
foreach([[1,2,8], [1,4,6], [2,5,6], [1,3,5]], v,
  my(a = R[v[1]], b = R[v[2]], c = R[v[3]], e = [-a, a, -b, b, -c, c]);
  my(r1 = ellrank(mkE1(a, b, c)), r2 = ellrank(mkE2(a, b, c)), t0 = getwalltime(), g = g2rank(e));
  printf("%s: E1 rank [%d,%d], E2 rank [%d,%d]  => J rank in [%d,%d];  descent bound %d (Sel dim %d, |S| = %d, %d s)\n",
         v, r1[1], r1[2], r2[1], r2[2], r1[1] + r2[1], r1[2] + r2[2], g[1], g[2], #g[3], (getwalltime() - t0) \ 1000);
  printf("    local dims: %s\n", g[4]));
}
