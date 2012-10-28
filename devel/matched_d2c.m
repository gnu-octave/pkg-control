numP = [1];
denP = conv ([1, 1, 1], [1, 4, 6, 4, 1]);
P = tf (numP, denP);

Pc = d2c (c2d (P, 0.1, 'matched'), 'matched')