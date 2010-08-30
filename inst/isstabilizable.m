## Copyright (C) 1998, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007
##               Kai P. Mueller.
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
## @deftypefn {Function File} {@var{retval} =} isstabilizable (@var{sys})
## @deftypefnx {Function File} {@var{retval} =} isstabilizable (@var{sys}, @var{tol})
## @deftypefnx {Function File} {@var{retval} =} isstabilizable (@var{a}, @var{b})
## @deftypefnx {Function File} {@var{retval} =} isstabilizable (@var{a}, @var{b}, @var{tol})
## @deftypefnx {Function File} {@var{retval} =} isstabilizable (@var{a}, @var{b}, @var{[]}, @var{dflg})
## @deftypefnx {Function File} {@var{retval} =} isstabilizable (@var{a}, @var{b}, @var{tol}, @var{dflg})
## Logical check for system stabilizability (i.e., all unstable modes are controllable).
##
## @strong{Inputs}
## @table @var
## @item sys
## LTI system.
## @item a
## State transition matrix.
## @item b
## Input matrix.
## @item tol
## Optional tolerance for stability. Default value is 0.
## @item dflg = 0
## Matrices (a, b) are part of a continuous-time system. Default Value.
## @item dflg = 1
## Matrices (a, b) are part of a discrete-time system.
## @end table
##
## @strong{Outputs}
## @table @var
## @item retval = 0
## System is not stabilizable.
## @item retval = 1
## System is stabilizable.
## @end table
##
## @example
## @group
## Method
## - Calculate staircase form (SLICOT AB01OD)
## - Extract unobservable part of state transition matrix
## - Calculate eigenvalues of unobservable part
## - Check whether
##   real (ev) < -tol*(1 + abs (ev))   continuous-time
##   abs (ev) < 1 - tol                discrete-time
## @end group
## @end example
## @seealso{isdetectable, isstable, isctrb, isobsv}
## @end deftypefn

## Author: A. S. Hodel <a.s.hodel@eng.auburn.edu>
## Created: August 1993
## Updated by A. S. Hodel (scotte@eng.auburn.edu) Aubust, 1995 to use krylovb
## Updated by John Ingram (ingraje@eng.auburn.edu) July, 1996 to accept systems

## Adapted-By: Lukas Reichlin <lukas.reichlin@gmail.com>
## Date: October 2009
## Version: 0.2

function retval = isstabilizable (a, b = [], tol = [], dflg = 0)

  if (nargin < 1 || nargin > 4)
    print_usage ();
  elseif (isa (a, "lti"))  # system passed
    if (nargin > 2)
      print_usage ();
    endif
    tol = b;
    dflg = isdt (a);
    [a, b] = ssdata (a);
  elseif (nargin == 1)  # a,b arguments sent directly
    print_usage ();
  elseif (! issquare (a) || rows (a) != rows (b))
    error ("isstabilizable: a must be square and conformal to b")
  endif

  if (isempty (tol))
    tol = 0;
  endif

  ## controllability staircase form
  [ac, bc, u, ncont] = slab01od (a, b, tol);

  ## extract uncontrollable part of staircase form
  uncont_idx = ncont+1 : rows (a);
  auncont = ac(uncont_idx, uncont_idx);

  ## calculate poles of uncontrollable part
  eigw = eig (auncont);

  ## check whether uncontrollable poles are stable
  if (dflg)
    retval = all (abs (eigw) < 1 - tol);
  else
    retval = all (real (eigw) < -tol*(1 + abs (eigw)));
  endif

endfunction
