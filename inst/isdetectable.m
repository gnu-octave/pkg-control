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
## @deftypefn {Function File} {@var{bool} =} isdetectable (@var{sys})
## @deftypefnx {Function File} {@var{bool} =} isdetectable (@var{sys}, @var{tol})
## @deftypefnx {Function File} {@var{bool} =} isdetectable (@var{a}, @var{c})
## @deftypefnx {Function File} {@var{bool} =} isdetectable (@var{a}, @var{c}, @var{tol})
## @deftypefnx {Function File} {@var{bool} =} isdetectable (@var{a}, @var{c}, @var{[]}, @var{dflg})
## @deftypefnx {Function File} {@var{bool} =} isdetectable (@var{a}, @var{c}, @var{tol}, @var{dflg})
## Logical test for system detectability. All unstable modes must be observable or
## all unobservable states must be stable. Uses SLICOT AB01OD by courtesy of NICONET e.V.
## <http://www.slicot.org>
##
## @strong{Inputs}
## @table @var
## @item sys
## LTI system.
## @item a
## State transition matrix.
## @item c
## Measurement matrix.
## @item tol
## Optional tolerance for stability. Default value is 0.
## @item dflg = 0
## Matrices (@var{a}, @var{c}) are part of a continuous-time system. Default Value.
## @item dflg = 1
## Matrices (@var{a}, @var{c}) are part of a discrete-time system.
## @end table
##
## @strong{Outputs}
## @table @var
## @item bool = 0
## System is not detectable.
## @item bool = 1
## System is detectable.
## @end table
##
## See @command{isstabilizable} for description of computational method.
## @seealso{isstabilizable, isstable, isctrb, isobsv}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.2

function bool = isdetectable (a, c = [], tol = [], dflg = 0)

  if (nargin < 1 || nargin > 4)
    print_usage ();
  elseif (isa (a, "lti"))  # isdetectable (sys), isdetectable (sys, tol)
    if (nargin > 2)
      print_usage ();
    endif
    tol = c;
    dflg = isdt (a);
    [a, b, c] = ssdata (a);
  elseif (nargin == 1)  # isdetectable (a, c, ...)
    print_usage ();
  endif

  bool = isstabilizable (a.', c.', tol, dflg);  # arguments checked inside

endfunction

