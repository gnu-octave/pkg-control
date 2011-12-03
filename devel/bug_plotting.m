s = tf ('s');
G = 1/(s^2+1)

[mag, pha] = bode (G)  % this works OK

figure
bode (G)               % fltk and gnuplot both crash