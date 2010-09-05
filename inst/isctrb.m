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
## @deftypefn {Function File} {[@var{bool}, @var{u}] =} isctrb (@var{sys})
## @deftypefnx {Function File} {[@var{bool}, @var{u}] =} isctrb (@var{sys}, @var{tol})
## @deftypefnx {Function File} {[@var{bool}, @var{u}] =} isctrb (@var{a}, @var{b})
## @deftypefnx {Function File} {[@var{bool}, @var{u}] =} isctrb (@var{a}, @var{b}, @var{tol})
## Logical check for system controllability.
## Uses SLICOT AB01OD by courtesy of NICONET e.V.
## <http://www.slicot.org>
##
## @strong{Inputs}
## @table @var
## @item sys
## LTI model.
## @item a
## State transition matrix.
## @item b
## Input matrix.
## @item tol
## Optional roundoff parameter. Default value is zero.
## @end table
##
## @strong{Outputs}
## @table @var
## @item bool = 0
## System is not controllable.
## @item bool = 1
## System is controllable.
## @item u
## An orthogonal basis of the controllable subspace.
## @end table
##
## @seealso{isobsv}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.2

function [bool, u] = isctrb (a, b = 0, tol = 0)

  if (nargin < 1 || nargin > 3)
    print_usage ();
  elseif (isa (a, "lti"))  # isctrb (sys), isctrb (sys, tol)
    if (nargin > 2)
      print_usage ();
    endif
    tol = b;
    [a, b] = ssdata (a);
  elseif (nargin < 2)  # isctrb (a, b), isctrb (a, b, tol)
    print_usage ();
  elseif (isempty (a) || isempty (b) || rows (a) != rows (b) || ! issquare (a))
    error ("isctrb: a(%dx%d), b(%dx%d) not conformal",
            rows (a), columns (a), rows (b), columns (b));
  endif

  ## check tol dimensions
  if (! isscalar (tol))
    error ("isctrb: tol(%dx%d) must be a scalar",
            rows (tol), columns (tol));
  endif

  [ac, bc, u, ncont] = slab01od (a, b, tol);

  u = u(:, 1:ncont);

  bool = (ncont == rows (a));

endfunction
