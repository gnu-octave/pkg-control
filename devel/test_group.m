sys = WestlandLynx;

sys.outgroup.attitude = [2, 3];
sys.outgroup.rate = [6, 5, 4];

sys(2:5,:)

% [sys; sys]

% [sys, sys]

sys ("attitude",:)

P = ss (toeplitz (1:5));

P = mktito (P, 2, 2)

P({'Y2', 'Y1'}, 'U2')
P({'Y2', 'Y1'}, {'U2', 'U1', 'U2'})