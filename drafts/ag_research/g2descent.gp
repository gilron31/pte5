\\ 2-descent for Jacobians of genus-2 curves y^2 = prod_{i=1..6} (t - e_i) with all
\\ e_i rational (full rational 2-torsion), giving a PROVEN upper bound on the rank.
\\
\\ 1. Odd model: t = e6 + 1/s gives w^2 = prod_{i<6} ((e6 - e_i) s + 1); then
\\    X = c s (c = prod (e6 - e_i)), W = c^2 w gives W^2 = prod_{i<6} (X - a_i), a_i = -c/(e6 - e_i).
\\ 2. x - T map:  mu(P - oo) = (X - a_1, ..., X - a_5) in (Q*/Q*^2)^5, with the
\\    usual fix at Weierstrass points; the image lies in {prod = square}.
\\ 3. Local images at p in S (2, primes dividing disc) and at oo: spanned by mu of
\\    random local points until the known size is reached: #J(Q_p)/2J(Q_p) = 2^4
\\    (p odd), 2^6 (p = 2); J(R)/2J(R) = 2^2 (5 real roots).
\\ 4. Sel = { delta in <-1, S>^5, prod square : res_v(delta) in image_v for all v };
\\    rank J(Q) <= dim Sel - 4.
\\ Entry point: g2rank(e) for a vector e of 6 distinct rationals.

\\ ---- square classes as F_2 vectors ----
sqc(x, p) =
{
  my(v, u);
  if (p == 0, return([if (x < 0, 1, 0)]));
  v = valuation(x, p); u = x / p^v;
  if (p == 2,
    u = lift(Mod(numerator(u), 8) / Mod(denominator(u), 8));
    \\ u = (-1)^a 5^b mod 8
    return([v % 2, if (u == 3 || u == 7, 1, 0), if (u == 3 || u == 5, 1, 0)]));
  [v % 2, if (kronecker(lift(Mod(numerator(u), p) / Mod(denominator(u), p)), p) == -1, 1, 0)];
}
issq_p(x, p) = { if (x == 0, return(1)); my(c = sqc(x, p)); c == vector(#c); }

odd_model(e) =
{
  my(e6 = e[6], c = prod(i = 1, 5, e6 - e[i]));
  vector(5, i, -c / (e6 - e[i]));
}

\\ mu of an affine non-Weierstrass point with X-coordinate X, as a list of 5 rationals
mu_pt(a, X) =
{
  my(j = 0);
  for (i = 1, 5, if (X == a[i], j = i));
  if (!j, return(vector(5, i, X - a[i])));
  my(v = vector(5, i, X - a[i]));
  v[j] = prod(i = 1, 5, if (i == j, 1, a[j] - a[i]));
  v;
}
locvec(v, p) = concat(vector(#v, i, sqc(v[i], p)));

\\ local image at p (0 = oo): F_2 basis (columns) of the span of mu over Q_p-points
local_image(a, p, target) =
{
  my(gens = List(), M, r = 0, tries = 0, f = prod(i = 1, 5, 'X - a[i]));
  for (j = 1, 5, listput(gens, locvec(mu_pt(a, a[j]), p)));      \\ 2-torsion
  M = matimage(Mod(Mat(vector(#gens, j, gens[j]~)), 2)); r = #M;
  while (r < target && tries < 200000,
    tries++;
    my(X);
    if (p == 0, X = (random(2001) - 1000) / (1 + random(50)),
      my(k = random(7) - 3, A = random(p^6) - p^6 \ 2);
      X = A * p^k + if (random(3) == 0, a[1 + random(5)], 0));
    my(fx = subst(f, 'X, X));
    if (fx == 0, next);
    if (p == 0, if (fx < 0, next), if (!issq_p(fx, p), next));
    my(w = Mod(locvec(mu_pt(a, X), p)~, 2), M2 = if (r, matconcat([M, w]), Mat(w)));
    if (matrank(M2) > r, M = M2; r++));
  [M, r, tries];
}

\\ main: returns [rank upper bound, dim Sel, S, local dims]
g2rank(e) =
{
  my(a = odd_model(e), S, gl, nb, H, rows = List(), locdims = List());
  \\ bad primes: 2 and the primes of the root differences e_i - e_j (every prime of
  \\ bad reduction of the odd model divides one of them); factor each separately
  S = [2];
  for (i = 1, 6, for (j = i + 1, 6,
    my(d = e[i] - e[j]); S = setunion(S, factor(abs(numerator(d)) * denominator(d))[,1]~)));
  for (i = 1, 6, if (denominator(e[i]) > 1, S = setunion(S, factor(denominator(e[i]))[,1]~)));
  S = select(q -> q > 1, S);
  gl = concat([-1], S);               \\ generators of <-1, S> mod squares
  nb = #gl;
  \\ basis of H(S): delta_1..delta_4 free, delta_5 = product (so prod is a square)
  H = List();
  for (i = 1, 4, for (k = 1, nb,
    my(d = vector(5, j, 1)); d[i] = gl[k]; d[5] = gl[k]; listput(H, d)));
  foreach(concat([0], S), p,
    my(target = if (p == 0, 2, if (p == 2, 6, 4)));
    my(LI = local_image(a, p, target));
    listput(locdims, [p, LI[2], target]);
    if (LI[2] < target, print("WARNING: local image at ", p, " has dim ", LI[2], " < ", target));
    \\ annihilator of the local image, applied to the localization of each basis vector
    \\ W: n_v x r (columns span the image); Ann: rows y with y*W = 0; A: n_v x #H
    my(W = LI[1], Ann = matker(Mod(W~, 2))~, A);
    A = Mod(Mat(vector(#H, k, locvec(H[k], p)~)), 2);
    if (#Ann, my(C = Ann * A); for (i = 1, matsize(C)[1], listput(rows, C[i,]))));
  my(Mc = matrix(#rows, #H, i, j, rows[i][j]), ker = matker(Mc));
  [#ker - 4, #ker, S, Vec(locdims)];
}
