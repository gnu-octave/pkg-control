## Copyright (C) 2010   Lukas F. Reichlin
##
## This file is part of LTI Syncope.
##
## LTI Syncope is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## LTI Syncope is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn{Function File} {@var{x} =} lyap (@var{a}, @var{b})
## @deftypefnx{Function File} {@var{x} =} lyap (@var{a}, @var{b}, @var{c})
## @deftypefnx{Function File} {@var{x} =} lyap (@var{a}, @var{b}, @var{[]}, @var{e})
## Solve continuous-time Lyapunov or Sylvester equations.
## Uses SLICOT SB03MD, SB04MD and SG03AD by courtesy of NICONET e.V.
## <http://www.slicot.org>
##
## @example
## @group
## AX + XA' + B = 0      (Lyapunov Equation)
##
## AX + XB + C = 0       (Sylvester Equation)
##
## AXE' + EXA' + B = 0   (Generalized Lyapunov Equation)
## @end group
## @end example
##
## @seealso{lyapchol, dlyap, dlyapchol}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: January 2010
## Version: 0.2

function [x, scale] = lyap (a, b, c, e)

  scale = 1;

  switch (nargin)
    case 2  # Lyapunov equation
  
      if (! is_real_square_matrix (a, b))
        error ("lyap: a, b must be real and square");
      endif
  
      if (rows (a) != rows (b))
        error ("lyap: a, b must have the same number of rows");
      endif

      [x, scale] = slsb03md (a, -b, false);  # AX + XA' = -B

      ## x /= scale;  # 0 < scale <= 1
    
    case 3  # Sylvester equation
    
      if (! is_real_square_matrix (a, b))
        error ("lyap: a, b must be real and square");
      endif

      if (! is_real_matrix (c) || rows (c) != rows (a) || columns (c) != columns (b))
        error ("lyap: c must be a real (%dx%d) matrix", rows (a), columns (b));
      endif

      x = slsb04md (a, b, -c);  # AX + XB = -C

    case 4  # generalized Lyapunov equation
    
      if (! isempty (c))
        print_usage ();
      endif
      
      if (! is_real_square_matrix (a, b, e))
        error ("lyap: a, b, e must be real and square");
      endif
      
      if (rows (b) != rows (a) || rows (e) != rows (a))
        error ("lyap: a, b, e must have the same number of rows");
      endif
      
      if (! issymmetric (b))
        error ("lyap: b must be symmetric");
      endif

      [x, scale] = slsg03ad (a, e, -b, false);  # AXE' + EXA' = -B
      
      ## x /= scale;  # 0 < scale <= 1

    otherwise
      print_usage ();

  endswitch

  if (scale < 1)
    warning ("lyap: solution scaled by %g to prevent overflow", scale);
  endif

endfunction


## Lyapunov
%!shared X, X_exp
%! A = [1, 2; -3, -4];
%! Q = [3, 1; 1, 1];
%! X = lyap (A, Q);
%! X_exp = [ 6.1667, -3.8333;
%!          -3.8333,  3.0000];
%!assert (X, X_exp, 1e-4);

## Sylvester
%!shared X, X_exp
%! A = [2.0   1.0   3.0
%!      0.0   2.0   1.0
%!      6.0   1.0   2.0];
%!
%! B = [2.0   1.0
%!      1.0   6.0];
%!
%! C = [2.0   1.0
%!      1.0   4.0
%!      0.0   5.0];
%!
%! X = lyap (A, B, -C);
%!
%! X_exp = [-2.7685   0.5498
%!          -1.0531   0.6865
%!           4.5257  -0.4389];
%!
%!assert (X, X_exp, 1e-4);

## Generalized Lyapunov
%!shared X, X_exp
%! A = [  3.0     1.0     1.0
%!        1.0     3.0     0.0
%!        1.0     0.0     2.0];
%!
%! E = [  1.0     3.0     0.0
%!        3.0     2.0     1.0
%!        1.0     0.0     1.0];
%!
%! B = [-64.0   -73.0   -28.0
%!      -73.0   -70.0   -25.0
%!      -28.0   -25.0   -18.0];
%!
%! X = lyap (A.', -B, [], E.');
%!
%! X_exp = [-2.0000  -1.0000   0.0000
%!          -1.0000  -3.0000  -1.0000
%!           0.0000  -1.0000  -3.0000];
%!
%!assert (X, X_exp, 1e-4);

