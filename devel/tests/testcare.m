%{
    READ ( NIN, FMT = * ) N, JOB, DICO, HINV, TRANA, UPLO, SCAL, SORT,
     $                      FACT, LYAPUN

Program Data
 SB02RD EXAMPLE PROGRAM DATA
   2     A     C     D     N     U     N     S     N     O
   0.0   1.0
   0.0   0.0
   1.0   0.0
   0.0   2.0
   0.0   0.0
   0.0   1.0
Program Results
 SB02RD EXAMPLE PROGRAM RESULTS

 The solution matrix X is 
   2.0000   1.0000
   1.0000   2.0000

 Estimated separation =   0.4000

 Estimated reciprocal condition number =   0.1333

 Estimated error bound =   0.0000

%}

%A,Q,G

A = [   0.0   1.0
        0.0   0.0 ];
   
Q = [   1.0   0.0
        0.0   2.0 ];
   
G = [   0.0   0.0
        0.0   1.0 ];


R = eye (2);

[B, p] = chol (G, "lower");

X = care (A, B, Q, R)
   
