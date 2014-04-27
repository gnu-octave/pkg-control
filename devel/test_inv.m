clear all, close all, clc;

P = Boeing707
Pi = inv (P)

p1 = pole (P)
z1 = zero (P)

pi1 = pole (Pi)
zi1 = zero (Pi)

figure (1)
sigma (P, Pi)

% until here, everything looks fine.  p1 == zi1, z1 == pi1

G = tf (P)
Gi = inv (G)
GPi = tf (Pi)

p2 = pole (G)
z2 = zero (G)

% FIXME: wrong/additional poles/zeros after second conversion to state-space
%        (needed for computation of MIMO TF zeros/poles)
pi2 = pole (Gi)
zi2 = zero (Gi)

ppi2 = pole (GPi)
zpi2 = zero (GPi)

figure (2)
sigma (G, Gi)


% FIXME: too many states (causes additional poles/zeros)
H = ss (Gi)
Hp = ss (GPi)


