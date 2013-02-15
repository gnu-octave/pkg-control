%{
s = tf ("s");

s + 2

s - 2

[s, 2]

[s; 2]


a = tf ([1, -2], [1, 3]);
b = tf ([1, -4], [1, 5]);

a*b

feedback (a, b)
feedback (a, b, "+")
%}

num = {[1, 5, 7], [1]; [1, 7], [1, 5, 5]};
den = {[1, 5, 6], [1, 2]; [1, 8, 6], [1, 3, 2]};
c = tf (num, den);

d = tf (Boeing707);

%c*d

%c+c
%{
[c, c]

[c; c]
%}

numcel = {[1 2],[ 1]; [0], [1 -1]};
dencel = {[1 1], [1 2 1]; [0], [4 1]};
mm = tf (numcel, dencel)
mms = ss (mm) # <- this crashed