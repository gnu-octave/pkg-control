## Copyright (C) 2017-2018 Fabian Alexander Wilms <f.alexander.wilms@gmail.com>
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
## @deftypefn {Function File} {@var{k} =} acker (@var{A}, @var{b}, @var{p})
## Calculates the state feedback matrix of a completely controllable SISO system
## using Ackermann's formula
##
## Given the state-space system
## @tex
## $$ \dot x = Ax + bu $$
## @end tex
## @ifnottex
## @example
## @group
##       .
##       x = Ax + bu
## @end group
## @end example
## @end ifnottex
## 
## and the desired closed-loop behaviour
##
## @tex
## $$ p(s) = s^n + \alpha_1 s^{n-1} + ... + \alpha_n $$
## @end tex
## @ifnottex
## @example
## @group
##       p(s) = s^n + a1*s^(n-1) + ... + an
## @end group
## @end example
## @end ifnottex
##
## via the vector p containing the roots of the polynom, the state feedback
## matrix k is calculated:
##
## @tex
## $$ k = (k_1 k_2 ... k_n) $$
## @end tex
## @ifnottex
## @example
## @group
##       k = (k1 k2 ... kn)
## @end group
## @end example
## @end ifnottex
##
## @seealso{place}
## @end deftypefn

## This function uses equation (4) from the paper
## "Sliding mode control design based on Ackermann's formula",
## DOI: 10.1109/9.661072

## Author: Fabian Alexander Wilms <f.alexander.wilms@gmail.com>
## Created: May 2017
## Version: 0.1

function K = acker(A, b, p)

  if (nargin != 3)
    print_usage ();
  endif
  
  if (! is_real_square_matrix (A) || ! is_real_vector (b) || rows (A) != rows (b))
      error ("acker: matrix A and vector b not conformal");
  endif
  
  if (! isnumeric (p) || ! isvector (p) || isempty (p))  # p could be complex
    error ("acker: p must be a vector");
  endif

  p = poly (p);

  p_A = zeros (size(A)(1), 1);

  for i=1:size (A)(1)+1
    p_A = p_A + p(i) * A ^ (size (A)(2) - (i - 1));
  end

  var = zeros (1, size (ctrb (A,b)));
  var(size(ctrb(A,b))) = 1;

  K = var * inv (ctrb (A, b)) * p_A;

endfunction

%!test
%! # https://en.wikipedia.org/wiki/Ackermann's_formula#Example
%! A = [1 1; 1 2];
%! B = [1; 0];
%! P = roots ([1 11 30]);
%! K = acker(A,B,P);
%! assert (K, [14 57]);
