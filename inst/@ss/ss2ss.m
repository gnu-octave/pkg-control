## Copyright (C) 2017 Fabian Alexander Wilms <f.alexander.wilms@gmail.com>
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {@var{SYS_TRANSFORMED} =} ss2ss (@var{SYS}, @var{T})
## @deftypefn {Function File} {[@var{A_TRANSFORMED} @var{B_TRANSFORMED} @var{C_TRANSFORMED} @var{D_TRANSFORMED}] =} ss2ss (@var{A},  @var{B}, @var{C}, @var{D}, @var{T})
## Applies the similarity transformation T to a state-space model
##
## Given the state space model
## @tex
## $$ \dot x = Ax + Bu $$
## $$ y = Cx + Du $$
## @end tex
## @ifnottex
## @example
## @group
##       .
##       x = Ax + Bu
##       y = Cx + Du
## @end group
## @end example
## @end ifnottex
##
## and the transformation matrix T, which maps the state vector x to another
## coordinate system
## @tex
## $$ \bar{x} = Tx $$
## @end tex
## @ifnottex
## @example
## @group
##       _
##       x = Tx
## @end group
## @end example
## @end ifnottex
##
## the state-space model is transformed in a way that results in an equivalent
## state-space model which is based on the new state vector
##
## @tex
## $$ \dot \bar{x} = TAT^{-1}\bar{x} + TBu $$
## $$ y = CT^{-1}\bar{x} + Du $$
## @end tex
## @ifnottex
## @example
## @group
##       .
##       _            _ 
##       x = T A T^-1 x + T B u
##       .          _
##       y = C T^-1 x + D u
## @end group
## @end example
## @end ifnottex
## 
## Please note: In the literature, T may be defined inversely:
##
## @tex
## $$ \bar{x} = T^{-1}x $$
## @end tex
## @ifnottex
## @example
## @group
##       _    -1
##       x = T  x
## @end group
## @end example
## @end ifnottex
## References:
##
## Control System Design, page 484 by Goodwin, Graebe, Salgado, 2000
##
## https://de.mathworks.com/help/control/ref/ss2ss.html
##
## Attention: T as defined by Matlab is the inverse of T as defined by 
## Goodwin, Graebe, Salgado

## Author: Fabian Alexander Wilms <f.alexander.wilms@gmail.com>

function [first_out, second_out, third_out, fourth_out] = ss2ss(first_in, second_in, third_in, fourth_in, fifth_in)
  
  if (nargin != 2 && nargin != 5)
    print_usage ();
  endif
  
  switch nargin
    case 2
      [A,B,C,D] = ssdata(first_in);
      % Attention: T as defined by Matlab is the inverse of T as defined by 
      % Goodwin, Graebe, Salgado
      T = inv(second_in);
    case 5;
      A = first_in;
      B = second_in;
      C = third_in;
      D = fourth_in;
      % see above
      T = inv(fifth_in);
  endswitch
   
  A_TRANSFORMED = inv(T)*A*T;
  B_TRANSFORMED = inv(T)*B;
  C_TRANSFORMED = C*T;
  D_TRANSFORMED = D;
  
  % make number of output variables depend on nargin
  switch nargin
    case 2
      if nargout > 1
        error('Too many output arguments')  
      endif
      first_out = ss(A_TRANSFORMED,B_TRANSFORMED,C_TRANSFORMED,D_TRANSFORMED);
    case 5
      if nargout > 4
        error('Too many output arguments')  
      endif
      first_out = A_TRANSFORMED;
      second_out = B_TRANSFORMED;
      third_out = C_TRANSFORMED;
      fourth_out = D_TRANSFORMED;
  endswitch
  
endfunction

%!test
%! A = [1 2 3; 4 5 6; 7 8 9]
%! B = [1; 2; 3]
%! C = [-1 0 1]
%! D = 0
%! 
%! [T E] = eig(A)
%! 
%! original_system = ss(A,B,C,D)
%! transformed_system = ss2ss(original_system,T)
%! assert(transformed_system.a, [0.5755    2.006  -0.8179; 6.684    13.62   -2.585; -4.476   -7.826   0.8073],e-14)
%! assert(transformed_system.b, [-0.5789; -3.148; 1.631],e-14)
%! assert(transformed_system.c, [0.8912  -0.2231    1.112],e-14)
%! assert(transformed_system.d, 0)
%! retransformed_system = ss2ss(transformed_system,inv(T))
%! assert(original_system.a,retransformed_system.a,e-14)
%! assert(original_system.b,retransformed_system.b,e-14)
%! assert(original_system.c,retransformed_system.c,e-14)
%! assert(original_system.d,retransformed_system.d,e-14)