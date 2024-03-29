## Copyright (C) 2009-2016   Lukas F. Reichlin
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
##
## @strong{Equations}
## @example
## @group
## A U' U A'  -  U' U  +  B B'  =  0           (Lyapunov Equation)
##
## A U' U A'  -  E U' U E'  +  B B'  =  0      (Generalized Lyapunov Equation)
## @end group
## @end example
##
## @strong{Algorithm}@*
## Uses @uref{https://github.com/SLICOT/SLICOT-Reference, SLICOT SB03OD and SG03BD},
## Copyright (c) 2020, SLICOT, available under the BSD 3-Clause
## (@uref{https://github.com/SLICOT/SLICOT-Reference/blob/main/LICENSE,  License and Disclaimer}).
##
## @seealso{dlyap, lyap, lyapchol}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: January 2010
## Version: 0.2.1

function [u, scale] = dlyapchol (a, b, e)

  switch (nargin)
    case 2
      
      if (! is_real_square_matrix (a))
        ## error ("dlyapchol: a must be real and square");
        error ("dlyapchol: %s must be real and square", ...
                inputname (1));
      endif

      if (! is_real_matrix (b))
        ## error ("dlyapchol: b must be real")
        error ("dlyapchol: %s must be real", ...
                inputname (2))
      endif
  
      if (rows (a) != rows (b))
        ## error ("dlyapchol: a and b must have the same number of rows");
        error ("dlyapchol: %s and %s must have the same number of rows", ...
                inputname (1), inputname (2));
      endif

      [u, scale] = __sl_sb03od__ (a.', b.', true);

      ## NOTE: TRANS = 'T' not suitable because we need U' U, not U U'

    case 3

      if (! is_real_square_matrix (a, e))
        ## error ("dlyapchol: a, e must be real and square");
        error ("dlyapchol: %s, %s must be real and square", ...
                inputname (1), inputname (3));
      endif

      if (! is_real_matrix (b))
        ## error ("dlyapchol: b must be real");
        error ("dlyapchol: %s must be real", ...
                inputname (2));
      endif

      if (rows (b) != rows (a) || rows (e) != rows (a))
        ## error ("dlyapchol: a, b, e must have the same number of rows");
        error ("dlyapchol: %s, %s, %s must have the same number of rows", ...
                inputname (1), inputname (2), inputname (3));
      endif

      [u, scale] = __sl_sg03bd__ (a.', e.', b.', true);

      ## NOTE: TRANS = 'T' not suitable because we need U' U, not U U'

    otherwise
      print_usage ();

  endswitch

  if (scale < 1)
    warning ("dlyapchol: solution scaled by %g to prevent overflow\n", scale);
  endif

endfunction


%!shared U, U_exp
%!
%! A = [ 0.5000    0.5000   -0.5000
%!       0.0000    1.0000   -0.5000
%!       1.0000   -0.5000    0.0000 ];
%!
%! B = [ -2.0 ; 1.0 ; -1.0 ];
%!
%! U = dlyapchol (A, B);
%!
%! U_exp = [ 5.6042    5.1284    1.9562
%!           0.0000    3.4313   -1.6716
%!           0.0000    0.0000    2.5621 ];
%!
%!assert (U, U_exp, 1e-4);

%!shared U, U_exp
%!
%! A = [ 0.5000    0.5000   -0.5000
%!       0.0000    1.0000   -0.5000
%!       1.0000   -0.5000    0.0000 ];
%!
%! B = [ -2.0 ; 1.0 ; -1.0 ];
%!
%! E = [ 2.0000    0.0000    0.0000
%!       0.0000    2.0000    0.0000
%!       0.0000    0.0000    2.0000 ];
%!
%! U = dlyapchol (A, B, E);
%!
%! U_exp = [ 1.0488   -0.3546    0.4644
%!           0.0000    0.6560   -0.5187
%!           0.0000    0.0000    0.4099 ];
%!
%!assert (U, U_exp, 1e-4);
