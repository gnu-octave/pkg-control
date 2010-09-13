## Copyright (C) 2009 - 2010   Lukas F. Reichlin
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
## @deftypefn {Function File} {[@var{x}, @var{l}, @var{g}] =} dare (@var{a}, @var{b}, @var{q}, @var{r})
## @deftypefnx {Function File} {[@var{x}, @var{l}, @var{g}] =} dare (@var{a}, @var{b}, @var{q}, @var{r}, @var{s})
## Solve discrete-time algebraic Riccati equation (ARE).
## Uses SLICOT SB02OD by courtesy of NICONET e.V.
## <http://www.slicot.org>
##
## @strong{Inputs}
## @table @var
## @item a
## Real matrix (n-by-n).
## @item b
## Real matrix (n-by-m).
## @item q
## Real matrix (n-by-n).
## @item r
## Real matrix (m-by-m).
## @item s
## Optional real matrix (n-by-m). If @var{s} is not specified, a zero matrix is assumed.
## @end table
##
## @strong{Outputs}
## @table @var
## @item x
## Unique stabilizing solution of the discrete-time Riccati equation (n-by-n).
## @item l
## Closed-loop poles (n-by-1).
## @item g
## Corresponding gain matrix (m-by-n).
## @end table
##
## @example
## @group
##                           -1
## A'XA - X - A'XB (B'XB + R)   B'XA + Q = 0
##
##                                 -1
## A'XA - X - (A'XB + S) (B'XB + R)   (A'XB + S)' + Q = 0
##
##               -1
## G = (B'XB + R)   B'XA
##
##               -1
## G = (B'XB + R)   (B'XA + S')
##
## L = eig (A - B*G)
## @end group
## @end example
## @seealso{care, lqr, dlqr, kalman}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.4

function [x, l, g] = dare (a, b, q, r, s = [])

  ## TODO: Add SLICOT SG02AD (Solution of continuous- or discrete-time
  ##       algebraic Riccati equations for descriptor systems)

  ## TODO: extract feedback matrix g from SB02OD (and SG02AD)

  if (nargin < 4 || nargin > 5)
    print_usage ();
  endif

  if (! isreal (a) || ! issquare (a))
    error ("dare: a must be real and square");
  endif
  
  if (! isreal (b) || rows (a) != rows (b))
    error ("dare: b must be real and conformal to a");
  endif

  if (! isreal (q) || ! issquare (q))
    error ("dare: q must be real and square");
  endif

  if (! isreal (r) || ! issquare (r))
    error ("dare: r must be real and square");
  endif
  
  if (columns (r) != columns (b))
    error ("dare: (b, r) not conformable");
  endif

  if (! isempty (s) && (! isreal (s) || any (size (s) != size (b))))
    error ("dare: s(%dx%d) must be real and identically dimensioned with b(%dx%d)",
            rows (s), columns (s), rows (b), columns (b));
  endif
  
  ## check stabilizability
  if (! isstabilizable (a, b, [], 1))
    error ("dare: (a, b) not stabilizable");
  endif

  ## check positive semi-definiteness
  if (isempty (s))
    t = zeros (size (b));
  else
    t = s;
  endif

  m = [q, t; t.', r];

  if (isdefinite (m) < 0)
    error ("dare: require [q, s; s.', r] >= 0");
  endif

  ## solve the riccati equation
  if (isempty (s))
    [x, l] = slsb02od (a, b, q, r, b, true, false);
    
    g = (r + b.'*x*b) \ (b.'*x*a);  # gain matrix
  else
    [x, l] = slsb02od (a, b, q, r, s, true, true);

    g = (r + b.'*x*b) \ (b.'*x*a + s.');  # gain matrix
  endif

endfunction


%!shared x, l, g, xe, le, ge
%! a = [ 0.4   1.7
%!       0.9   3.8];
%!
%! b = [ 0.8
%!       2.1];
%!
%! c = [ 1  -1];
%!
%! r = 3;
%!
%! [x, l, g] = dare (a, b, c.'*c, r);
%!
%! xe = [ 1.5354    1.2623
%!        1.2623   10.5596];
%!
%! le = [-0.0022
%!        0.2454];
%!
%! ge = [ 0.4092    1.7283];
%!
%!assert (x, xe, 1e-4);
%!assert (sort (l), sort (le), 1e-4);
%!assert (g, ge, 1e-4);