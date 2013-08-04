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


sys = ss (A, B, C, D)

z1 = zero (sys, 'system')       % {-1, 2, -4} system zeros

z2 = zero (sys, 'invariant')    % {2, -1} invariant zeros

z3 = zero (sys, 'transmission') % {2} transmission zeros

z4 = zero (sys, 'output')       % {-1} output decoupling zeros

z5 = zero (sys, 'input')        % {-4} input decoupling zeros







%z3 = zero (minreal (sys2))
%z4 = zero (sminreal (sys2))


