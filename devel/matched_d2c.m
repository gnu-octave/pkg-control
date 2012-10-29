numP = [1];
denP = conv ([1, 1, 1], [1, 4, 6, 4, 1]);
P = tf (numP, denP);


Pd = c2d (P, 0.1, 'matched')
[z, k] = zero (Pd)
p = pole (Pd)


Pc = d2c (c2d (P, 0.1, 'matched'), 'matched')