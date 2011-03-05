P = ss (-2, 3, 4, 0)

Pi = inv (P)

p = pole (P)
[z, k] = zero (P)

pi = pole (Pi)
[zi, ki] = zero (Pi)


%{
P = dss (-2, 3, 4, 5, 1)
[z, k] = zero (P)
[z, k] = zero (ss (-2, 3, 4, 5))
%}