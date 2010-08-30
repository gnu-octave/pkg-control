## Copyright (C) 1993, 1994, 1995, 2000, 2002, 2003, 2004, 2005, 2006,
##               2007 Auburn University.  All rights reserved.
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
## @deftypefn {Function File} {@var{retval} =} isdetectable (@var{sys})
## @deftypefnx {Function File} {@var{retval} =} isdetectable (@var{sys}, @var{tol})
## @deftypefnx {Function File} {@var{retval} =} isdetectable (@var{a}, @var{c})
## @deftypefnx {Function File} {@var{retval} =} isdetectable (@var{a}, @var{c}, @var{tol})
## @deftypefnx {Function File} {@var{retval} =} isdetectable (@var{a}, @var{c}, @var{[]}, @var{dflg})
## @deftypefnx {Function File} {@var{retval} =} isdetectable (@var{a}, @var{c}, @var{tol}, @var{dflg})
## Logical test for system detectability (observability of unstable modes).
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
## Matrices (a, c) are part of a continuous-time system. Default Value.
## @item dflg = 1
## Matrices (a, c) are part of a discrete-time system.
## @end table
##
## @strong{Outputs}
## @table @var
## @item retval = 0
## System is not detectable.
## @item retval = 1
## System is detectable.
## @end table
##
## See @command{isstabilizable} for detailed description of
## arguments and computational method.
## @seealso{isstabilizable, isstable, isctrb, isobsv}
## @end deftypefn

## Author: A. S. Hodel <a.s.hodel@eng.auburn.edu>
## Created: August 1993
## Updated by John Ingram (ingraje@eng.auburn.edu) July 1996.

## Adapted-By: Lukas Reichlin <lukas.reichlin@gmail.com>
## Date: October 2009
## Version: 0.2

function retval = isdetectable (a, c = [], tol = [], dflg = 0)

  if (nargin < 1 || nargin > 4)
    print_usage ();
  elseif (isa (a, "lti"))  # system passed
    if (nargin > 2)
      print_usage ();
    endif
    tol = c;
    dflg = isdt (a);
    [a, b, c] = ssdata (a);
  elseif (nargin == 1)  # a,c arguments sent directly
    print_usage ();     # a,c dimension checked inside isstabilizable
  endif

  retval = isstabilizable (a.', c.', tol, dflg);

endfunction

