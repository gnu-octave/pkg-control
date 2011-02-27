P = ss (-2, 3, 4, 5);

u = [zeros(1, 20), ones(1, 20)];
t = 5;

figure (1)
lsim (P, u, t)
% lsim (P, u)

figure (2)
lsim (c2d (P, 0.1), u)

figure (3)
lsim (c2d (P, 0.1), u, 3.9)

% figure (4)
% lsim (c2d (P, 0.1), u, 0:0.1:3.9)
% lsim (Boeing707, ones (10, 2), 0:9)


t = 0:0.1:3.9;
figure (4)
lsim (P, u, t)

figure (5)
lsim (P, u, t+100)  # This one works, ...

figure (6)
lsim (P, u, t+30)   # ... but I have no idea why this one fails!