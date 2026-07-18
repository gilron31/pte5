# Quadric parameterization

## Usecase 1

$$
x^2 + y^2 = z^2 \\

P = (0, 1, 1) \\

l = (0, 1, 1) + t(1, a ,b)
$$

Solving:

$$

t^2 + (1+at)^2 = (1+bt)^2 \\

t^2 + 1 + 2at + a^2t^2 - 1 - 2bt - b^2t^2 = 0 \\

t(1 + a^2 - b^2) + 2(a-b) = 0 \\

t = \frac{2(b-a)}{(1 + a^2 - b^2)} \\

(x, y, z) = (\frac{2(b-a)}{(1 + a^2 - b^2)}, 1  + a\frac{2(b-a)}{(1 + a^2 - b^2)}, 1 + b\frac{2(b-a)}{(1 + a^2 - b^2)}) \\

(x, y, z) = (\frac{2(b-a)}{(1 + a^2 - b^2)},  \frac{2ab-a^2 + 1 - b^2}{(1 + a^2 - b^2)}, \frac{b^2-2ab + 1 + a^2 }{(1 + a^2 - b^2)}) \\

(x, y, z) = (\frac{2(b-a)}{(1 + a^2 - b^2)},  \frac{-(a-b)^2 +  1 }{(1 + a^2 - b^2)}, \frac{(a-b)^2 + 1  }{(1 + a^2 - b^2)}) \\

(x, y, z) = (2(b-a),  -(a-b)^2 +  1 , (a-b)^2 + 1  ) \\

$$

is it true?

$$

(4d^2 + d^4 - 2d^2 + 1 - d^4 - 2d^2 - 1 = 0)

$$

Again with fuller parameterization

$$
x^2 + y^2 = z^2 \\

P = (0, 1, 1) \\

l = (0, 1, 1) + t(a, b ,c)
$$

Solving:

$$

a^2t^2 + (1+bt)^2 = (1+ct)^2 \\

a^2t^2 + 1 + 2bt + b^2t^2 - 1 - 2ct - c^2t^2 = 0 \\

t(a^2 + b^2 - c^2) + 2(b-c) = 0 \\

t = \frac{2(c-b)}{(a^2 + b^2 - c^2)} \\

(x, y, z) = (at, 1  + bt, 1 + ct) \\

(x, y, z) = (a2(c-b), a^2 + b^2 - c^2  + b2(c-b), a^2 + b^2 - c^2 + c2(c-b)) \\

(x, y, z) = (a2(c-b), a^2 - b^2 - c^2  + 2bc, a^2 + b^2 + c^2 - 2bc) \\

(x, y, z) = (a2(c-b), a^2 - (b-c)^2, a^2 + (b-c)^2) \\

d = c-b \\

(x, y, z) = (2ad, a^2 - d^2, a^2 + d^2) \\

$$

## Usecase 2

$$
x^2 + y^2 = z^2 + w^2 \\

P = (0, 1, 0, 1) \\

l = P + t(a, b, c, d) \\

a^2t^2 + (1 + bt)^2 - (ct)^2 - (1 + dt)^2 = 0 \\

t^2(a^2 + b^2 - c^2 - d^2) + 2t(b - d) = 0 \\

t = \frac{2(d - b)}{(a^2 + b^2 - c^2 - d^2)} \\

(x, y, z, w) = (at, 1+bt, ct, 1+dt) \\

m = d - b \\

(x, y, z, w) = (2am, a^2 - c^2 - m^2 , 2cm, a^2 - c^2 + m^2) \\

$$
