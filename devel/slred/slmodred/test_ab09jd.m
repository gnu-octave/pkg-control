% AB09JD EXAMPLE PROGRAM DATA (Continuous system)
%  6     1     1     2   0     0   0.0  1.E-1  1.E-14    V   N   I   C    S   A

a =  [ -3.8637   -7.4641   -9.1416   -7.4641   -3.8637   -1.0000
        1.0000,         0         0         0         0         0
             0    1.0000         0         0         0         0
             0         0    1.0000         0         0         0
             0         0         0    1.0000         0         0
             0         0         0         0    1.0000         0 ];

b =  [       1
             0
             0
             0
             0
             0 ];

c =  [       0         0         0         0         0         1 ];

d =  [       0 ];

av = [  0.2000   -1.0000
        1.0000         0 ];

bv = [       1
             0 ];

cv = [ -1.8000         0 ];

dv = [       1 ];


[ar, br, cr, dr] = slab09jd (a, b, c, d, av, bv, cv, dv, [], [], [], [], 0, 0.0, \
                             1, 0, 2, 0, 0, 1, \
                             1e-1, 1e-14)

%{
 The reduced state dynamics matrix Ar is 
  -0.2391   0.3072   1.1630   1.1967
  -2.9709  -0.2391   2.6270   3.1027
   0.0000   0.0000  -0.5137  -1.2842
   0.0000   0.0000   0.1519  -0.5137

 The reduced input/state matrix Br is 
  -1.0497
  -3.7052
   0.8223
   0.7435

 The reduced state/output matrix Cr is 
  -0.4466   0.0143  -0.4780  -0.2013

 The reduced input/output matrix Dr is 
   0.0219
%}
