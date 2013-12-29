A = [1 0; 0 1];
B = [1; 0];
C = [1 0];
D = 0;
E = [1 0; 0 0];


sys = dss (A, B, C, D, E)
t = 0:0.1:1;
u = sin (t);
lsim (sys, u, t)
