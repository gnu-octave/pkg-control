sys1 = frd (1, 1)

sys2 = frd (2+5i, 3, 0.1)

sys3 = frd (ones (2, 3, 5), 1:5)

sys4 = frd

sys5 = frd ([], [])

sys6 = frd ([1, 2; 3, 4], logspace (-1, 1, 10))



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


%sys1 = frd (reshape (1:60, 3, 4, []), 1:5)
%sys1.'
P_frd.';

%sys2 = frd (reshape (1:45, 3, 3, []), 1:5)
%inv (sys2)
inv (P_frd);


P_ss = Boeing707;
T_ss = feedback (P_ss);

P_frd = frd (P_ss);
T_frd = feedback (P_frd);

figure (5)
subplot (2, 1, 1)
% sigma (T_ss)
sigma (T_ss, T_frd.w)
subplot (2, 1, 2)
sigma (T_frd)

figure (6)
subplot (2, 1, 1)
% sigma (T_ss + P_ss)
sigma (T_ss + P_ss, T_frd.w)
subplot (2, 1, 2)
sigma (T_frd + P_frd)

figure (7)
subplot (2, 1, 1)
% sigma (T_ss * P_ss)
sigma (T_ss * P_ss, T_frd.w)
subplot (2, 1, 2)
sigma (T_frd * P_frd)


k1 = frd (ss ([1, 2; 3, 4]));
% k2 = frd (tf ([1, 2; 3, 4]))
k2 = frd (tf (5));

frd (ss (5)) + 4;


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