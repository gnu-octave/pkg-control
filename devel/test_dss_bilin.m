sys = inv (Boeing707);

d = c2d (sys, 2, "bilin")

c = d2c (d, "bilin")

d.e+d.a
rcond (d.e + d.a)
