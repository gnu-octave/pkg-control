a = [ -1.0   0.0   0.0
      -2.0  -2.0  -2.0
      -1.0   0.0  -3.0 ];
  
b = [  1.0   0.0   0.0
       0.0   2.0   1.0 ].';
   
   
c = [  0.0   2.0   1.0
       1.0   0.0   0.0 ];


[ac, bc, cc, z, ncont] = sltb01ud (a, b, c, 0.0)

a = [ 1     1
      4    -2 ];

b = [ 1    -1
      1    -1 ];

c = [ 1     0
      0     1 ];

[ac, bc, cc, z, ncont] = sltb01ud (a, b, c, 0.0)

%{
 The transformed state dynamics matrix of a controllable realization is 
  -3.0000   2.2361
   0.0000  -1.0000

 and the dimensions of its diagonal blocks are 
  2

 The transformed input/state matrix B of a controllable realization is 
   0.0000  -2.2361
   1.0000   0.0000

 The transformed output/state matrix C of a controllable realization is 
  -2.2361   0.0000
   0.0000   1.0000

 The controllability index of the transformed system representation =  1

 The similarity transformation matrix Z is 
   0.0000   1.0000   0.0000
  -0.8944   0.0000  -0.4472
  -0.4472   0.0000   0.8944
%}
%{
A =
     1     1
     4    -2

B =
     1    -1
     1    -1

C =
     1     0
     0     1
and locate the uncontrollable mode.

[Abar,Bbar,Cbar,T,k]=ctrbf(A,B,C)

Abar =
   -3.0000         0
   -3.0000    2.0000

Bbar =
    0.0000    0.0000
    1.4142   -1.4142

Cbar =
   -0.7071    0.7071
    0.7071    0.7071

T =
   -0.7071    0.7071
    0.7071    0.7071
k =
     1     0
%}