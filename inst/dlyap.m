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
## @deftypefn{Function File} {@var{x} =} dlyap (@var{a}, @var{b})
## @deftypefnx{Function File} {@var{x} =} dlyap (@var{a}, @var{b}, @var{c})
## @deftypefnx{Function File} {@var{x} =} dlyap (@var{a}, @var{b}, @var{[]}, @var{e})
## Solve discrete-time Lyapunov or Sylvester equations.
##
## @strong{Equations}
## @example
## @group
## AXA' - X + B = 0      (Lyapunov Equation)
##
## AXB - X + C = 0       (Sylvester Equation)
##
## AXA' - EXE' + B = 0   (Generalized Lyapunov Equation)
## @end group
## @end example
##
## @strong{Algorithm}@*
## Uses @uref{https://github.com/SLICOT/SLICOT-Reference, SLICOT SB03MD, SB04QD and SG03AD},
## Copyright (c) 2020, SLICOT, available under the BSD 3-Clause
## (@uref{https://github.com/SLICOT/SLICOT-Reference/blob/main/LICENSE,  License and Disclaimer}).
##
## @seealso{dlyapchol, lyap, lyapchol}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: January 2010
## Version: 0.2.1

function [x, scale] = dlyap (a, b, c, e)

  scale = 1;

  switch (nargin)
    case 2                                     # Lyapunov equation

      if (! is_real_square_matrix (a, b))
        ## error ("dlyap: a, b must be real and square");
        error ("dlyap: %s, %s must be real and square", ...
                inputname (1), inputname (2));
      endif
  
      if (rows (a) != rows (b))
        ## error ("dlyap: a, b must have the same number of rows");
        error ("dlyap: %s, %s must have the same number of rows", ...
                inputname (1), inputname (2));
      endif

      if issymmetric (b)

        ## The 'normal' case where b is symmetric
        [x, scale] = __sl_sb03md__ (a, -b, true);     # AXA' - X = -B
        ## x /= scale;                                # 0 < scale <= 1

      else

        ## b is non-symmetric, solve as Sylvester equation
        x = __sl_sb04qd__ (-a, a', b);    # AXB - X = -C  (A = a, B = a', C = b)

      endif
  
  
    case 3                                     # Sylvester equation
  
      if (! is_real_square_matrix (a, b))
        ## error ("dlyap: a, b must be real and square");
        error ("dlyap: %s, %s must be real and square", ...
                inputname (1), inputname (2));
      endif

      if (! is_real_matrix (c) || rows (c) != rows (a) || columns (c) != columns (b))
        ## error ("dlyap: c must be a real (%dx%d) matrix", rows (a), columns (b));
        error ("dlyap: %s must be a real (%dx%d) matrix", ...
                rows (a), columns (b), inputname (3));
      endif

      x = __sl_sb04qd__ (-a, b, c);                 # AXB' - X = -C

    case 4                                     # generalized Lyapunov equation
          
      if (! isempty (c))
        print_usage ();
      endif
      
      if (! is_real_square_matrix (a, b, e))
        ## error ("dlyap: a, b, e must be real and square");
        error ("dlyap: %s, %s, %s must be real and square", ...
                inputname (1), inputname (2), inputname (4));
      endif
      
      if (rows (b) != rows (a) || rows (e) != rows (a))
        ## error ("dlyap: a, b, e must have the same number of rows");
        error ("dlyap: %s, %s, %s must have the same number of rows", ...
                inputname (1), inputname (2), inputname (4));
      endif
      
      if (! issymmetric (b))
        ## error ("dlyap: b must be symmetric");
        error ("dlyap: %s must be symmetric", ...
                inputname (2));
      endif

      [x, scale] = __sl_sg03ad__ (a, e, -b, true);  # AXA' - EXE' = -B
      
      ## x /= scale;                           # 0 < scale <= 1

    otherwise
      print_usage ();

  endswitch

  if (scale < 1)
    warning ("dlyap: solution scaled by %g to prevent overflow\n", scale);
  endif

endfunction


## Lyapunov
%!shared X, X_exp
%! A = [3.0   1.0   1.0
%!      1.0   3.0   0.0
%!      0.0   0.0   3.0];
%!
%! B = [25.0  24.0  15.0
%!      24.0  32.0   8.0
%!      15.0   8.0  40.0];
%!
%! X = dlyap (A.', -B);
%!
%! X_exp = [2.0000   1.0000   1.0000
%!          1.0000   3.0000   0.0000
%!          1.0000   0.0000   4.0000];
%!
%!assert (X, X_exp, 1e-4);

## Lyapunov with non-symmetric B
%!shared X, X_exp
%! A = [3.0   1.0   1.0
%!      1.0   3.0   0.0
%!      0.0   0.0   3.0];
%!
%! B = [1.0   0.2   2.0
%!      0.5   0.5   1.0
%!      1.0  -2.0   1.0];
%!
%! X = dlyap (A.', -B);
%!
%! X_exp = [0.1390   -0.0514    0.1831
%!         -0.0086    0.0676    0.0422
%!          0.2005   -0.3233   -0.0362];
%!
%!assert (X, X_exp, 1e-4);

## Sylvester
%!shared X, X_exp
%! A = [1.0   2.0   3.0 
%!      6.0   7.0   8.0 
%!      9.0   2.0   3.0];
%!
%! B = [7.0   2.0   3.0 
%!      2.0   1.0   2.0 
%!      3.0   4.0   1.0];
%!
%! C = [271.0   135.0   147.0
%!      923.0   494.0   482.0
%!      578.0   383.0   287.0];
%!
%! X = dlyap (-A, B, C);
%!
%! X_exp = [2.0000   3.0000   6.0000
%!          4.0000   7.0000   1.0000
%!          5.0000   3.0000   2.0000];
%!
%!assert (X, X_exp, 1e-4);

## Generalized Lyapunov
%!shared X, X_exp
%! A = [3.0   1.0   1.0
%!      1.0   3.0   0.0
%!      1.0   0.0   2.0];
%!
%! E = [1.0   3.0   0.0
%!      3.0   2.0   1.0
%!      1.0   0.0   1.0];
%!
%! B = [ -3.0  -10.0   7.0
%!      -10.0  -14.0  -2.0
%!        7.0   -2.0   9.0 ];
%!
%! X = dlyap (A, B, [], E);
%!
%! X_exp = [-2.0000   -1.0000   0.0000
%!          -1.0000   -3.0000  -1.0000
%!           0.0000   -1.0000  -3.0000];
%!
%!assert (X, X_exp, 1e-4);
