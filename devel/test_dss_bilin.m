sys = inv (Boeing707);

d = c2d (sys, 0.2, "bilin")

c = d2c (d, "bilin")
