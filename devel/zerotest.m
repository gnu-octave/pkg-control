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

P_tf = tf (P)
Pi_tf = tf (Pi)

p_tf = pole (P_tf)
[z_tf, k_tf] = zero (P_tf)

pi_tf = pole (Pi_tf)
[zi_tf, ki_tf] = zero (Pi_tf)