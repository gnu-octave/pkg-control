ss_sys = WestlandLynx;
tf_sys = tf (ss_sys);

w = logspace (-4, 3, 500);

figure (1)
sigma (ss_sys, w)

figure (2)
sigma (tf_sys, w)


sys = ss (0, 1, 1, 0)
tf (sys)
