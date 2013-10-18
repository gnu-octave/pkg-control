sys = WestlandLynx;

sys.outgroup.attitude = [2, 3];
sys.outgroup.rate = [6, 5, 4];

sys(2:5,:)

% [sys; sys]

% [sys, sys]

sys ("attitude",:)

P = ss (toeplitz (1:5));

P = mktito (P, 2, 2)

P({'V', 'Z'}, 'U')
P({'V', 'Z'}, {'U', 'W', 'U'})

P({"z3", "V"},:)


xperm (sys, {'q', 'p'})
