## Examples from <http://ch.mathworks.com/help/control/ref/lsim.html>
clear all, close all, clc;

figure (1)
H = [tf([2 5 1],[1 2 3]);tf([1 -1],[1 1 5])];
[u,t] = gensig('square',4,10,0.1);
lsim(H,u,t)


figure (2)
w2 = 62.83^2;
h = tf(w2,[1 2 w2]);
t = 0:0.1:5;                % vector of time samples
u = (rem(t,1) >= 0.5);        % square wave values
lsim(h,u,t)


figure (3)
dt = 0.016;
ts = 0:dt:5;
us = (rem(ts,1) >= 0.5);
hd = c2d(h,dt);
lsim(hd,us,ts)

