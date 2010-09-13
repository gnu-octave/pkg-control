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
## @deftypefn {Function File} {[@var{bool}, @var{u}] =} isobsv (@var{sys})
## @deftypefnx {Function File} {[@var{bool}, @var{u}] =} isobsv (@var{sys}, @var{tol})
## @deftypefnx {Function File} {[@var{bool}, @var{u}] =} isobsv (@var{a}, @var{c})
## @deftypefnx {Function File} {[@var{bool}, @var{u}] =} isobsv (@var{a}, @var{c}, @var{tol})
## Logical check for system observability.
## Uses SLICOT AB01OD by courtesy of NICONET e.V.
## <http://www.slicot.org>
##
## @strong{Inputs}
## @table @var
## @item sys
## LTI model.
## @item a
## State transition matrix.
## @item c
## Measurement matrix.
## @item tol
## Optional roundoff parameter. Default value is zero.
## @end table
##
## @strong{Outputs}
## @table @var
## @item bool = 0
## System is not observable.
## @item bool = 1
## System is observable.
## @item u
## An orthogonal basis of the observable subspace.
## @end table
##
## @seealso{isctrb}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.2.1

function [bool, u] = isobsv (a, c = [], tol = [])

  if (nargin < 1 || nargin > 3)
    print_usage ();
  elseif (isa (a, "lti"))  # isobsv (sys), isobsv (sys, tol)
    if (nargin > 2)
      print_usage ();
    endif
    tol = c;
    [a, b, c] = ssdata (a);
  elseif (nargin < 2)  # isobsv (a, c), isobsv (a, c, tol)
    print_usage ();
  endif

  [bool, u] = isctrb (a.', c.', tol);

endfunction

