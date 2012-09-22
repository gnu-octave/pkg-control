load tfs.dat

figure (1)
step2 (T_AH, T_opt)

figure (2)
step2 (T_AH, '-.r', T_opt, '-b')

figure (3)
step2 (T_AH, c2d (T_opt, 1))
