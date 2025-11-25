mm1 = tf(-.50405, [1 .98752] , -1)

figure (1)
step (mm1)

figure (2)
step (mm1, 450)



mm2 = tf(-.50405, [1 .98752] , 1)

figure (3)
step (mm2)

[y, t, x] = step (mm2);

sy = size (y)
st = size (t)
sx = size (x)