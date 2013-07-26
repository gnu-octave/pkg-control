A = diag ([1, 1, 3, -4, -1, 3]);

B = [   0,  -1
       -1,   0
        1,  -1
        0,   0
        0,   1
       -1,  -1  ];
       
C = [   1,  0,  0,  1,  0,  0
        0,  1,  0,  1,  0,  1
        0,  0,  1,  0,  0,  1   ];
        
D = zeros (3, 2);


sys1 = ss (A, B, C, D)

z1 = zero (sys1)    % {2, -1} invariant zeros

z2 = tzero (sys1)   % {2} transmission zeros

z3 = szero (sys1)   % {-1, 2, -4} system zeros

sys4 = ss (A, B, zeros(0,columns(A)), zeros(0, columns(B)));
z4 = zero (sys4)    % {-4} input decoupling zeros

sys5 = ss (A, zeros(rows(A), 0), C, zeros(rows(C),0));
z5 = zero (sys5)    % {-1} output decoupling zeros





%z3 = zero (minreal (sys2))
%z4 = zero (sminreal (sys2))


