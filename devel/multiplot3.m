load tfs.dat

figure (1)
step (T_AH, T_opt)

figure (2)
step (T_AH, '-.r', T_opt, '-b')

figure (3)
step (T_AH, c2d (T_opt, 1))

figure (4)
lsim (T_AH, T_opt, ones(3000,1), 30)