## Copyright (C) 2009   Lukas F. Reichlin
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
## @deftypefn {Function File} {[@var{g}, @var{x}, @var{l}] =} dlqr (@var{sys}, @var{q}, @var{r})
## @deftypefnx {Function File} {[@var{g}, @var{x}, @var{l}] =} dlqr (@var{sys}, @var{q}, @var{r}, @var{s})
## @deftypefnx {Function File} {[@var{g}, @var{x}, @var{l}] =} dlqr (@var{a}, @var{b}, @var{q}, @var{r})
## @deftypefnx {Function File} {[@var{g}, @var{x}, @var{l}] =} dlqr (@var{a}, @var{b}, @var{q}, @var{r}, @var{s})
## Return linear-quadratic state-feedback gain matrix g for a LTI system as well as
## the solution x of the associated riccati equation and the closed-loop poles l.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: November 2009
## Version: 0.1

function [g, x, l] = dlqr (a, b, q, r = [], s = [])

  if (nargin < 3 || nargin > 5)
    print_usage ();
  endif

  if (isa (a, "lti"))
    s = r;
    r = q;
    q = b;
    [a, b, c, d, tsam] = ssdata (a);
  elseif (nargin < 4)
    print_usage ();
  else
    tsam = 1;  # any value > 0 could be used here
  endif

  if (tsam > 0)
    [x, l, g] = dare (a, b, q, r, s);
  else
    [x, l, g] = care (a, b, q, r, s);
  endif

endfunction