sysc = Boeing707
sysd = d2c (c2d (sysc, 1))


sys = ss (1, 2, 3, 4, 5)
d2c (sys)


sysd = c2d (sysc, 2, "prewarp", 100)
syse = d2c (sysd, "prewarp", 100)
sysc


sysc_tf = tf (sysc)
sysd_tf = c2d (sysc_tf, 0.2)
syse_tf = d2c (sysd_tf)
