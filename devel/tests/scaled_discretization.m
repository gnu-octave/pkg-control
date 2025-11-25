P = BMWengine ("unscaled")



% Ps = prescale (P)

Ps = P;

P.scaled = true;    % damits bei c2d nicht skaliert wird


x0 = 2:6;


figure (1)
initial (P, x0)

figure (2)
initial (Ps, x0)

figure (3)
initial (P, Ps, x0)







#{

A = [ -2.8    2.0   -1.8
      -2.4   -2.0    0.8
       1.1    1.7   -1.0 ];

B = [ -0.8    0.5    0
       0      0.7    2.3
      -0.3   -0.1    0.5 ];

C = [ -0.1    0     -0.3
       0.9    0.5    1.2
       0.1   -0.1    1.9 ];

D = [ -0.5    0      0
       0.1    0      0.3
      -0.8    0      0   ];

x_0 = [1, 2, 3];

sysc = ss (A, B, C, D);

sysc
prescale (sysc)

[yc, tc, xc] = initial (sysc, x_0, 0.2, 0.1);
initial_c = [yc, tc, xc];

sysd = c2d (sysc, 2);

[yd, td, xd] = initial (sysd, x_0, 4);
initial_d = [yd, td, xd];

## expected values computed by the "dark side"

yc_exp = [ -1.0000    5.5000    5.6000
           -0.9872    5.0898    5.7671
           -0.9536    4.6931    5.7598 ];

tc_exp = [  0.0000
            0.1000
            0.2000 ];

xc_exp = [  1.0000    2.0000    3.0000
            0.5937    1.6879    3.0929
            0.2390    1.5187    3.0988 ];

initial_c_exp = [yc_exp, tc_exp, xc_exp];

yd_exp = [ -1.0000    5.5000    5.6000
           -0.6550    3.1673    4.2228
           -0.5421    2.6186    3.4968 ];

td_exp = [  0
            2
            4 ];

xd_exp = [  1.0000    2.0000    3.0000
           -0.4247    1.5194    2.3249
           -0.3538    1.2540    1.9250 ];

initial_d_exp = [yd_exp, td_exp, xd_exp];

assert (initial_c, initial_c_exp, 1e-4)
assert (initial_d, initial_d_exp, 1e-4)
#}


%{
shared initial_c, initial_c_exp, initial_d, initial_d_exp

A = [ -2.8    2.0   -1.8
      -2.4   -2.0    0.8
       1.1    1.7   -1.0 ];

B = [ -0.8    0.5    0
       0      0.7    2.3
      -0.3   -0.1    0.5 ];

C = [ -0.1    0     -0.3
       0.9    0.5    1.2
       0.1   -0.1    1.9 ];

D = [ -0.5    0      0
       0.1    0      0.3
      -0.8    0      0   ];

x_0 = [1, 2, 3];

sysc = ss (A, B, C, D);

[yc, tc, xc] = initial (sysc, x_0, 0.2, 0.1);
initial_c = [yc, tc, xc];

sysd = c2d (sysc, 2);

[yd, td, xd] = initial (sysd, x_0, 4);
initial_d = [yd, td, xd];

## expected values computed by the "dark side"

yc_exp = [ -1.0000    5.5000    5.6000
           -0.9872    5.0898    5.7671
           -0.9536    4.6931    5.7598 ];

tc_exp = [  0.0000
            0.1000
            0.2000 ];

xc_exp = [  1.0000    2.0000    3.0000
            0.5937    1.6879    3.0929
            0.2390    1.5187    3.0988 ];

initial_c_exp = [yc_exp, tc_exp, xc_exp];

yd_exp = [ -1.0000    5.5000    5.6000
           -0.6550    3.1673    4.2228
           -0.5421    2.6186    3.4968 ];

td_exp = [  0
            2
            4 ];

xd_exp = [  1.0000    2.0000    3.0000
           -0.4247    1.5194    2.3249
           -0.3538    1.2540    1.9250 ];

initial_d_exp = [yd_exp, td_exp, xd_exp];

assert (initial_c, initial_c_exp, 1e-4)
assert (initial_d, initial_d_exp, 1e-4)
%}