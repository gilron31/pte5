\\ Ranks of elliptic quotients of the bottom-lift curve over Symes' E_4 (roots below).
\\ 4-subsets: y^2 = prod (t + e), e in S (every 91st of the 1820 subsets);
\\ triples: E1 = the u = t^2 quotient, E2 = the u*v quotient; ord(P_N) = 0 means infinite order.
\\ Run: gp -q lift_ranks.gp
read("lift_lib.gp");
R = [131, 343, 11, 367, 101, 353, 77, 359];
N = 2 * 67405;
signed = concat(R, -R);
{
print("--- 4-subsets of the 16 signed roots (sample of 20) ---");
cnt = 0;
forsubset([16, 4], v,
  cnt++;
  if (cnt % 91 == 1,
    my(al = vector(4, i, -signed[v[i]]), E = quartic_E(al), rk = ellrank(E));
    print(Vec(v), "  rank in [", rk[1], ",", rk[2], "]  tors ", elltors(E)[2])));
}
{
print("--- triples: E1 (u=t^2 quotient) and E2 ---");
forsubset([8, 3], v,
  my(a = R[v[1]], b = R[v[2]], c = R[v[3]], e1 = mkE1(a, b, c), e2 = mkE2(a, b, c));
  my(r1 = ellrank(e1), r2 = ellrank(e2), PN, ip);
  PN = [N, sqrtint((N - a^2) * (N - b^2) * (N - c^2))];
  ip = if (ellisoncurve(e1, PN), ellorder(e1, PN), "off");
  printf("%s  E1 rank [%d,%d] tors %s ord(P_N)=%s | E2 rank [%d,%d] tors %s\n",
         Vec(v), r1[1], r1[2], elltors(e1)[2], ip, r2[1], r2[2], elltors(e2)[2]));
}
