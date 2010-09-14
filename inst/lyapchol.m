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
## @deftypefn{Function File} {@var{u} =} lyapchol (@var{a}, @var{b})
## @deftypefnx{Function File} {@var{u} =} lyapchol (@var{a}, @var{b}, @var{e})
## Compute Cholesky factor of continuous-time Lyapunov equations.
## Uses SLICOT SB03OD and SG03BD by courtesy of NICONET e.V.
## <http://www.slicot.org>
##
## @example
## @group
## A U' U  +  U' U A'  +  B B'  =  0           (Lyapunov Equation)
##
## A U' U E'  +  E U' U A'  +  B B'  =  0      (Generalized Lyapunov Equation)
## @end group
## @end example
##
## @seealso{lyap, dlyap, dlyapchol}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: January 2010
## Version: 0.2

function [u, scale] = lyapchol (a, b, e)

  switch (nargin)
    case 2
      
      if (! isreal (a) || isempty (a) || ! issquare (a))
        error ("lyapchol: a must be real and square");
      endif

      if (! isreal (b) || isempty (b))
        error ("lyapchol: b must be real and square")
      endif
  
      if (rows (a) != rows (b))
        error ("lyapchol: a and b must have the same number of rows");
      endif

      [u, scale] = slsb03od (a.', b.', false);

      ## NOTE: TRANS = 'T' not suitable because we need U' U, not U U'

    case 3

      if (! isreal (a) || isempty (a) || ! issquare (a))
        error ("lyapchol: a must be real and square");
      endif
      
      if (! isreal (e) || isempty (e) || ! issquare (e))
        error ("lyapchol: e must be real and square");
      endif

      if (! isreal (b) || isempty (b))
        error ("lyapchol: b must be real");
      endif

      if (rows (b) != rows (a) || rows (e) != rows (a))
        error ("lyapchol: a, b, e not conformal");
      endif

      [u, scale] = slsg03bd (a.', e.', b.', false);

      ## NOTE: TRANS = 'T' not suitable because we need U' U, not U U'

    otherwise
      print_usage ();

  endswitch

  if (scale < 1)
    warning ("lyapchol: solution scaled by %g to prevent overflow", scale);
  endif

endfunction


%!shared U, U_exp, X, X_exp
%!
%! A = [ -1.0  37.0 -12.0 -12.0
%!       -1.0 -10.0   0.0   4.0
%!        2.0  -4.0   7.0  -6.0
%!        2.0   2.0   7.0  -9.0 ].';
%!
%! B = [  1.0   2.5   1.0   3.5
%!        0.0   1.0   0.0   1.0
%!       -1.0  -2.5  -1.0  -1.5
%!        1.0   2.5   4.0  -5.5
%!       -1.0  -2.5  -4.0   3.5 ].';
%!
%! U = lyapchol (A, B);
%!
%! X = U.' * U;  # use lyap at home!
%!
%! U_exp = [  1.0000   0.0000   0.0000   0.0000
%!            3.0000   1.0000   0.0000   0.0000
%!            2.0000  -1.0000   1.0000   0.0000
%!           -1.0000   1.0000  -2.0000   1.0000 ].';
%!
%! X_exp = [  1.0000   3.0000   2.0000  -1.0000
%!            3.0000  10.0000   5.0000  -2.0000
%!            2.0000   5.0000   6.0000  -5.0000
%!           -1.0000  -2.0000  -5.0000   7.0000 ];
%!
%!assert (U, U_exp, 1e-4);
%!assert (X, X_exp, 1e-4);

%!shared U, U_exp, X, X_exp
%!
%! A = [ -1.0    3.0   -4.0
%!        0.0    5.0   -2.0
%!       -4.0    4.0    1.0 ].';
%!
%! E = [  2.0    1.0    3.0
%!        2.0    0.0    1.0
%!        4.0    5.0    1.0 ].';
%!
%! B = [  2.0   -1.0    7.0 ].';
%!
%! U = lyapchol (A, B, E);
%!
%! U_exp = [  1.6003  -0.4418  -0.1523
%!            0.0000   0.6795  -0.2499
%!            0.0000   0.0000   0.2041 ];
%!
%!assert (U, U_exp, 1e-4);