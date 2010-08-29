## Copyright (C) 1993, 1994, 1995, 2000, 2002, 2004, 2005, 2006, 2007
##               Auburn University.  All rights reserved.
##
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{retval}, @var{u}] =} isctrb (@var{sys}, @var{tol})
## @deftypefnx {Function File} {[@var{retval}, @var{u}] =} isctrb (@var{a}, @var{b}, @var{tol})
## Logical check for system controllability.
## Uses SLICOT AB01OD by courtesy of NICONET e.V.
## <http://www.slicot.org>
##
## @strong{Inputs}
## @table @var
## @item sys
## system data structure
## @item a
## @itemx b
## @var{n} by @var{n}, @var{n} by @var{m} matrices, respectively
## @item tol
## optional roundoff parameter. Default value: 0
## @end table
##
## @strong{Outputs}
## @table @var
## @item retval
## Logical flag; returns true (1) if the system @var{sys} or the
## pair (@var{a}, @var{b}) is controllable, whichever was passed as input
## arguments.
## @item u
## @var{u} is an orthogonal basis of the controllable subspace.
## @end table
## @seealso{isobsv}
## @end deftypefn

## Author: A. S. Hodel <a.s.hodel@eng.auburn.edu>
## Created: August 1993
## Updated by A. S. Hodel (scotte@eng.auburn.edu) Aubust, 1995 to use krylovb
## Updated by John Ingram (ingraje@eng.auburn.edu) July, 1996 for packed systems

## Adapted-By: Lukas Reichlin <lukas.reichlin@gmail.com>
## Date: October 2009
## Version: 0.2

function [retval, U] = isctrb (A, B = 0, tol = 0)

  if (nargin < 1 || nargin > 3)
    print_usage ();
  elseif (isa (A, "lti"))  # isctrb (sys), isctrb (sys, tol)
    if (nargin > 2)
      print_usage ();
    endif
    tol = B;
    [A, B] = ssdata (A);
  elseif (nargin < 2)  # isctrb (A, B), isctrb (A, B, tol)
    print_usage ();
  endif

  ## check tol dimensions
  if (! isscalar (tol))
    error ("isctrb: tol(%dx%d) must be a scalar",
            rows (tol), columns (tol));
  endif

  ## check dimensions compatibility
  if (isempty (A) || isempty (B) || rows (A) != rows (B) || ! issquare (A))
    error ("isctrb: A(%dx%d), B(%dx%d)",
            rows (A), columns (A), rows (B), columns (B));
  endif

  [Ac, Bc, U, ncont] = slab01od (A, B, tol);
U = U(:, 1:ncont);
  retval = (ncont == rows (A));

endfunction
