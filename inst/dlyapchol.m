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
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with LTI Syncope.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn{Function File} {@var{u} =} dlyapchol (@var{a}, @var{b})
## @deftypefnx{Function File} {@var{u} =} dlyapchol (@var{a}, @var{b}, @var{e})
## Compute Cholesky factor of discrete-time Lyapunov equations.
## Uses SLICOT SB03OD and SG03BD by courtesy of NICONET e.V.
## <http://www.slicot.org>
##
## @example
## @group
## A U' U A'  -  U' U  +  B B'  =  0           (Lyapunov Equation)
##
## A U' U A'  -  E U' U E'  +  B B'  =  0      (Generalized Lyapunov Equation)
## @end group
## @end example
##
## @seealso{dlyap, lyap, lyapchol}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: January 2010
## Version: 0.2

function [u, scale] = dlyapchol (a, b, e)

  switch (nargin)
    case 2
      
      if (! is_real_square_matrix (a))
        error ("dlyapchol: a must be real and square");
      endif

      if (! is_real_matrix (b))
        error ("dlyapchol: b must be real")
      endif
  
      if (rows (a) != rows (b))
        error ("dlyapchol: a and b must have the same number of rows");
      endif

      [u, scale] = slsb03od (a.', b.', true);

      ## NOTE: TRANS = 'T' not suitable because we need U' U, not U U'

    case 3

      if (! is_real_square_matrix (a, e))
        error ("dlyapchol: a, e must be real and square");
      endif

      if (! is_real_matrix (b))
        error ("dlyapchol: b must be real");
      endif

      if (rows (b) != rows (a) || rows (e) != rows (a))
        error ("dlyapchol: a, b, e must have the same number of rows");
      endif

      [u, scale] = slsg03bd (a.', e.', b.', true);

      ## NOTE: TRANS = 'T' not suitable because we need U' U, not U U'

    otherwise
      print_usage ();

  endswitch

  if (scale < 1)
    warning ("dlyapchol: solution scaled by %g to prevent overflow", scale);
  endif

endfunction


## TODO: add tests
