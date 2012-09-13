load tfs.dat

figure (1)
bode2 (C_AH, C_opt)

figure (2)
bode2 (5*C_AH, frd (C_AH), frd (C_opt))

% bode2 (5*C_AH, frd (C_AH, 1:10), frd (C_opt, 11:20))

figure (3)
bode2 (5*C_AH, '*r', C_AH, ':b', C_opt, '-.k')