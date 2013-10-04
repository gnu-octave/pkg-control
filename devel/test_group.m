sys = WestlandLynx;

sys.outgroup.attitude = [2, 3];
sys.outgroup.rate = [6, 5, 4];

sys(2:5,:)

% [sys; sys]

% [sys, sys]

sys ("attitude",:)
