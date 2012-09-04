a = [ -1.0   0.0   0.0
       0.0  -2.0   0.0
       0.0   0.0  -3.0 ];

b = [  0.0   1.0  -1.0
       1.0   1.0   0.0 ].';

c = [  0.0   1.0   1.0
       1.0   1.0   1.0 ];

d = [  1.0   0.0
       0.0   1.0 ];

[p, m] = size (d);
md = 4

[gn, gd, ign, igd] = __sl_tb04bd__ (a, b, c, d)

num = reshape (gn, md, p, m)
den = reshape (gd, md, p, m)

numc = mat2cell (num, md, ones(1,p), ones(1,m))
denc = mat2cell (den, md, ones(1,p), ones(1,m))

numc = squeeze (numc)
denc = squeeze (denc)

ignc = mat2cell (ign, ones(1,p), ones(1,m));
igdc = mat2cell (igd, ones(1,p), ones(1,m));

num = cellfun (@(x, y) x(1:y+1), numc, ignc, "uniformoutput", false);
den = cellfun (@(x, y) x(1:y+1), denc, igdc, "uniformoutput", false);


tf (num, den)

%{
numc{1,2,1}
denc{1,2,1}
numc{1,1,2}
denc{1,1,2}
%}


%num(:, 1, 1)
%den(:, 1, 1)
%num(:, 1, 2)
%den(:, 1, 2)

%num = mat2cell (gn, p, m)
%den = mat2cell (gd, p, m)

%num = gn(1:ign+1)
%den = gd(1:igd+1)
%{
[gn, gd, ign, igd] = __sl_tb04bd__ (-2, 3, 4, 5)


% for i = 1 : size (gn, 1)
%   for j = 1 : size (gn, 2)
%     gn(i, j, :)
%     gd(i, j, :)
%   endfor
% endfor


P = tf (1, [1 5 11 14 11 5 1]);

S = ss (P);

[num, den, ign, igd] = __sl_tb04bd__ (S.a, S.b, S.c, S.d)

P



[num, den, ign, igd] = __sl_tb04bd__ (0, 1, 1, 0)


sys = WestlandLynx;

[num, den, ign, igd] = __sl_tb04bd__ (sys.a, sys.b, sys.c, sys.d);
%}
