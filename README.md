# Toy PLONK

This toy application construct proof that `3 ^ 2 + 4 ^ 2 == 5 ^ 2`. Using small prime group (mod 101) and small elliptic curve.

```
cargo test -- circuit::circuit --nocapture
```

Will print:

```
circuit:
f_a = 1 + 13 * x + 3 * x ^ 2 + 3 * x ^ 3
f_b = 7 + 3 * x + 14 * x ^ 2 + 13 * x ^ 3
f_c = 6 + 5 * x + 11 * x ^ 2 + 4 * x ^ 3
q_l = 13 + x + 4 * x ^ 2 + 16 * x ^ 3
q_r = 13 + x + 4 * x ^ 2 + 16 * x ^ 3
q_0 = 16
q_m = 5 + 16 * x + 13 * x ^ 2 + x ^ 3
q_c = 0
s_1 = 7 + 13 * x + 10 * x ^ 2 + 6 * x ^ 3
s_2 = 4 + 13 * x ^ 2 + x ^ 3
s_3 = 6 + 7 * x + 3 * x ^ 2 + 14 * x ^ 3

proof:
a(s) = { 91, 66 }
b(s) = { 26, 45 }
c(s) = { 91, 35 }
zet(s) = { 32, 59 }
t_lo(s) = { 12, 32 }
t_mid(s) = { 26, 45 }
t_hi(s) = { 91, 66 }
w_zeta(s) = { 91, 35 }
w_zeta_omega(s) = { 65, 98 }
a(zeta) = 15
b(zeta) = 13
c(zeta) = 5
sigma1(zeta) = 1
sigma2(zeta) = 12
r(zeta) = 15
zet(zeta * omega) = 15

zeta ^ n - 1 = 12
L1(zeta) = 5
t(zeta) = 1
d = inf
f = { 68, 27 }
e = { 1, 2 }
lhs = { 32, 42 }
rhs = { 12, 69 }
l = { 93, 76u }, r = { 93, 76u }
```
