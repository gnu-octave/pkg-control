sys = frd (1, 1)

sys = frd (2+5i, 3, 0.1)

sys = frd (ones (2, 3, 5), 1:5)

frd

frd ([], [])



P_ss = Boeing707;
w = logspace (-2, 2, 500);
H = freqresp (P_ss, w);
P_frd = frd (H, w);

figure (1)
subplot (2, 1, 1)
sigma (P_ss)
subplot (2, 1, 2)
sigma (P_frd)

figure (2)
subplot (2, 1, 1)
sigma (P_ss, [], 1)
subplot (2, 1, 2)
sigma (P_frd, [], 1)

figure (3)
subplot (2, 1, 1)
sigma (P_ss, [], 2)
subplot (2, 1, 2)
sigma (P_frd, [], 2)

figure (4)
subplot (2, 1, 1)
sigma (P_ss, [], 3)
subplot (2, 1, 2)
sigma (P_frd, [], 3)


sys1 = frd (reshape (1:60, 3, 4, []), 1:5)
sys1.'

sys2 = frd (reshape (1:45, 3, 3, []), 1:5)
inv (sys2)

%{
P = frd ((1:45)*(1-2i), logspace (-2,3,45))

figure (1)
sigma (P)
sigma (P, [], 0)

figure (2)
sigma (P, [], 1)

figure (3)
sigma (P, [], 2)

figure (4)
sigma (P, [], 3)
%}