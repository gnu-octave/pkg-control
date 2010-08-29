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
## Return unique stabilizing solution x of the discrete-time
## Riccati equation as well as the closed-loop poles l and the
## corresponding gain matrix g.
## Uses SLICOT SB02OD by courtesy of NICONET e.V.
## <http://www.slicot.org>
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
## @end group
## @end example
## @seealso{care, lqr, dlqr, kalman}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.2

function [x, l, g] = dare (a, b, q, r, s = [])

  ## TODO: Add SLICOT SG02AD (Solution of continuous- or discrete-time
  ##       algebraic Riccati equations for descriptor systems)

  ## TODO: Check stabilizability and controllability more elegantly
  ##       (without incorporating cross terms into a and q)

  if (nargin < 4 || nargin > 5)
    print_usage ();
  endif

  [brows, bcols] = size (b);

  if (! issquare (a))
    error ("dare: a is not square");
  endif

  if (! issquare (q))
    error ("dare: q is not square");
  endif

  if (! issquare (r))
    error ("dare: r is not square");
  endif
  
  if (rows (a) != brows)
    error ("dare: a, b are not conformable");
  endif
  
  if (columns (r) != bcols)
    error ("dare: b, r are not conformable");
  endif

  ## incorporate cross term into a and q
  if (isempty (s))
    ao = a;
    qo = q;
  else
    [srows, scols] = size (s);  % [n2, m2]

    if (srows != brows || scols != bcols)
      error ("dare: s (%dx%d) must be identically dimensioned with b (%dx%d)",
              srows, scols, brows, bcols);
    endif

    ao = a - (b/r)*s.';
    qo = q - (s/r)*s.';
  endif
  
  ## check stabilizability
  if (! isstabilizable (ao, b, [], 1))
    error ("dare: a and b not stabilizable");
  endif

  ## check detectability
  dflag = isdetectable (ao, qo, [], 1);

  if (dflag == 0)
    warning ("dare: (a,q) not detectable");
  elseif (dflag == -1)
    error ("dare: (a,q) has non-minimal modes near unit circle");
  endif

  ## to allow lqe design, don't force
  ## qo to be positive semi-definite

  ## checking positive definiteness
  if (isdefinite (r) <= 0)
    error ("dare: r must be positive definite");
  endif

  ## solve the riccati equation
  if (isempty (s))
    ## unique stabilizing solution
    x = slsb02od (a, b, q, r, b, true, false);
    
    ## corresponding gain matrix
    g = (r + b.'*x*b) \ (b.'*x*a);
  else
    ## unique stabilizing solution
    x = slsb02od (a, b, q, r, s, true, true);
    
    ## corresponding gain matrix
    g = (r + b.'*x*b) \ (b.'*x*a + s.');
  endif

  ## closed-loop poles
  l = eig (a - b*g);
  
  ## TODO: use alphar, alphai and beta from SB02OD

endfunction


## TODO: add a test