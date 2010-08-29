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
## @deftypefn {Function File} {@var{retval} =} isdetectable (@var{a}, @var{c}, @var{tol}, @var{dflg})
## @deftypefnx {Function File} {@var{retval} =} isdetectable (@var{sys}, @var{tol})
## Test for detectability (observability of unstable modes) of (@var{a}, @var{c}).
##
## Returns 1 if the system @var{a} or the pair (@var{a}, @var{c}) is
## detectable, 0 if not, and -1 if the system has unobservable modes at the
## imaginary axis (unit circle for discrete-time systems).
##
## @strong{See} @command{is_stabilizable} for detailed description of
## arguments and computational method.
## @seealso{isstabilizable, size, rows, columns, length, ismatrix, isscalar, isvector}
## @end deftypefn

## Author: A. S. Hodel <a.s.hodel@eng.auburn.edu>
## Created: August 1993
## Updated by John Ingram (ingraje@eng.auburn.edu) July 1996.

## Adapted-By: Lukas Reichlin <lukas.reichlin@gmail.com>
## Date: October 2009
## Version: 0.1

function retval = isdetectable (a, c = [], tol = [], dflg = 0)

  if (nargin < 1)
    print_usage ();
  elseif (isa (a, "lti"))
    ## system form
    if (nargin == 2)
      tol = c;
    elseif (nargin > 2)
      print_usage ();
    endif
    dflg = isdt (a);
    [a, b, c] = ssdata (a);
  else
    if (nargin > 4 || nargin == 1)
      print_usage ();
    endif
  endif

  retval = isstabilizable (a.', c.', tol, dflg);

endfunction

