## Copyright (C) 2009, 2010   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{sys} =} minreal (@var{sys})
## @deftypefnx {Function File} {@var{sys} =} minreal (@var{sys}, @var{tol})
## Minimal realization or zero-pole cancellation of LTI models.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.2

function sys = minreal (sys, tol = "def")

  if (nargin > 2)  # nargin == 0 not possible because minreal is inside @lti
    print_usage ();
  endif

  if (! is_real_scalar (tol) && tol != "def")
    error ("minreal: second argument must be a real scalar");
  endif

  sys = __minreal__ (sys, tol);

endfunction