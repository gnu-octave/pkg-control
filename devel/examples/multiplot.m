load tfs.dat

figure (1)
bode (C_AH, C_opt)

figure (2)
bode (5*C_AH, frd (C_AH), frd (C_opt))

% bode2 (5*C_AH, frd (C_AH, 1:10), frd (C_opt, 11:20))

figure (3)
bode (5*C_AH, '*r', C_AH, 'xb', C_opt, 'ok')


figure (4)
nyquist (C_AH, C_opt)


figure (5)
nyquist (C_AH, "xr", C_opt, "ob")
legend ("Test C_AH", "Test C_opt")