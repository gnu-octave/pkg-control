feedback (tf (1), tf (-1))

%{
Transfer function 'ans' from input 'u1' to output ...

      1
 y1:  -
      0

Static gain.
%}


% prevent this in @tf/__sys_connect__.m around line 59;
%     sys.den(:) = den{1,1} * den{2,2}  -  M(1,2) * M(2,1) * num{1,1} * num{2,2};
% den must not be zero!


sys = feedback (ss (1), ss (-1))

%{
sys.e =
       x1  x2
   x1   0   0
   x2   0   0

sys.a =
       x1  x2
   x1  -1   1
   x2   1  -1

sys.b =
       u1
   x1   1
   x2   0

sys.c =
       x1  x2
   y1   1   0

sys.d =
       u1
   y1   0

Static gain.
%}

% Because E is a zero matrix, sys can be simplified to
% y = [D - C A^-1 B] u

gain = sys.d - sys.c / sys.a * sys.b

% However, the A matrix was not invertible

rcond (sys.a)

% So we have an algebraic loop if
% E is zero and A is not invertible
% ! any (sys.e(:)) && rcond (sys.a) < eps 


% Now let's see this system:
G = dss (1, 2, 3, 4, 0)

% G.a is invertible and G.e is zero, so
Gmin = G.d - G.c / G.a * G.b


% Conclusions:
% - @ss/__set__.m should not accept all-zero E matrices if A is not invertible,
%   because the system contains algebraic loops.  Therefore it would not be
%   possible to enter such systems in the 'dss' function.
%   (For the same reason we refuse zero denominators in 'tf')
%
% - @ss/__minreal__.m (and possibly sminreal) should use the formula
%   y = [D - C A^-1 B] u  if E is zero and A is invertible
%
% - @ss/__sys_connect__.m should error out because of algebraic loops if
%   E is zero and A is not invertible
%   ! any (sys.e(:)) && rcond (sys.a) < eps
%   otherwise (zero E and invertible A),
%   use  y = [D - C A^-1 B] u  and set  sys.e = []

% I'll have to consult my pillow :-)


