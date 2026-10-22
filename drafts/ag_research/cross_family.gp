\\ Cross-family test: points Q of family 2's (s:t:u) curve whose phi-image lies on
\\ family 1's image curve F1(J1, J2) = 0; rational candidates, then preimages on family 1.
default(parisizemax, 2^32);
L1 = t^2 + u^2 + 6;
L2 = 2*t^4 - 4*t^2*u^2 + 16*t^2 + 2*u^4 + 16*u^2 + 32;
L3 = 16*t^6 - 16*t^4*u^2 - 128*t^4 - 16*t^2*u^4 + 256*t^2*u^2 + 256*t^2 + 16*u^6 - 128*u^4 + 256*u^2;
G2 = -25*t^6 - 12*t^5*u^2 + 70*t^5*u + 28*t^4*u^3 + 29*t^4*u^2 - 40*t^4*u - 400*t^4 + 28*t^3*u^4 - 68*t^3*u^3 - 120*t^3*u^2 - 200*t^3*u + 500*t^3 - 12*t^2*u^5 + 29*t^2*u^4 - 120*t^2*u^3 + 800*t^2*u^2 - 100*t^2*u + 500*t^2 + 70*t*u^5 - 40*t*u^4 - 200*t*u^3 - 100*t*u^2 - 3000*t*u - 25*u^6 - 400*u^4 + 500*u^3 + 500*u^2 + 10000;
G1 = 16*t^6 - 32*t^5*u - 92*t^4*u^2 + 61*t^4 + 236*t^3*u^3 - 436*t^3*u - 92*t^2*u^4 + 806*t^2*u^2 - 500*t^2 - 32*t*u^5 - 436*t*u^3 + 1500*t*u + 16*u^6 + 61*u^4 - 500*u^2;
F1 = 4918547515869140625*J1^12-59022570190429687500*J1^11-4281103757812500000*J1^10*J2+326892400754508984375*J1^10+42811037578125000000*J1^9*J2+1397352266550000000*J1^8*J2^2-1104763100562667968750*J1^9-183903103140597112500*J1^8*J2-11178818132400000000*J1^7*J2^2-202709235467520000*J1^6*J2^3+2536046767443831663125*J1^8+443759923249776900000*J1^7*J2+48769918226432150400*J1^6*J2^2+1216255412805120000*J1^5*J2^3+11027382409433088*J1^4*J2^4-4162032079562247055000*J1^7-655729615987865913120*J1^6*J2-136116055504992902400*J1^5*J2^2+1639466655925761792*J1^4*J2^3-44109529637732352*J1^3*J2^4+5001367129334382437868*J1^6+598640495680318878720*J1^5*J2+241301135729268707652*J1^4*J2^2-14666236042403847168*J1^3*J2^3+1113803780385300480*J1^2*J2^4-4427944605883623357208*J1^5-310699300092815690804*J1^4*J2-266435444688588814608*J1^3*J2^2+24691228254188093952*J1^2*J2^3-2139388501495136256*J1*J2^4+83029702847496192*J2^5+2862576089724425555264*J1^4+55566319831639038416*J1^3*J2+176784105976581283608*J1^2*J2^2-16806636656088173568*J1*J2^3+1024516883293507584*J2^4-1316041094272324866336*J1^3+28814533017685140880*J1^2*J2-64432743843699344400*J1*J2^2+4134246272875452672*J2^3+407900906299896867072*J1^2-17938545882084948672*J1*J2+9910088717681350148*J2^2-76442618460404367360*J1+2959378437280942080*J2+6544227548115763200;
\\ N(t,u) = L1^44 * F1(L2/L1^2, L3/L1^4), then strip powers of L1
N = substvec(F1, ['J1, 'J2], [L2/L1^2, L3/L1^4]) * L1^44;
N = numerator(N);
while (type(N / L1) == "t_POL" && denominator(N / L1) == 1, N = N / L1);
printf("pulled-back curve: total degree in u %d, in t %d\n", poldegree(N, 'u), poldegree(N, 't));
t0 = getwalltime();
R = polresultant(G2, N, 'u);
printf("resultant in t: degree %d (%d s)\n", poldegree(R, 't), (getwalltime() - t0) \ 1000);
t0 = getwalltime();
F = factor(R);
writebin("cross_family_F.bin", F);
printf("factored (%d s): %d distinct factors; degrees %s\n", (getwalltime() - t0) \ 1000, #F[,1], vecsort(vector(#F[,1], i, poldegree(F[i,1], 't))));
\\ rational t-roots and the rational points of family 2's curve above them
{
cand = List();
for (i = 1, #F[,1],
  if (poldegree(F[i,1], 't) != 1, next);
  my(tt = -polcoef(F[i,1], 0, 't) / polcoef(F[i,1], 1, 't));
  my(g = gcd(subst(G2, 't, tt), subst(N, 't, tt)));
  if (poldegree(g, 'u) < 1, next);
  my(G = factor(g)[,1]);
  for (k = 1, #G,
    if (poldegree(G[k], 'u) != 1, next);
    my(uu = -polcoef(G[k], 0, 'u) / polcoef(G[k], 1, 'u));
    listput(cand, [tt, uu])));
printf("rational points of family 2's curve with phi-image on family 1's image: %d\n", #cand);
print(cand);
}
