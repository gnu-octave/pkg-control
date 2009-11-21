## Copyright (C) 1997, 2000, 2002, 2004, 2005, 2006, 2007 Kai P. Mueller
## Copyright (C) 2009   Lukas F. Reichlin
## Copyright (C) 2009 Luca Favatella <slackydeb@gmail.com>
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
## @deftypefn {Function File} {@var{ob} =} obsv (@var{sys})
## @deftypefnx {Function File} {@var{ob} =} obsv (@var{a}, @var{c})
## Build observability matrix:
## @iftex
## @tex
## $$ O_b = \left[ \matrix{  C       \cr
##                           CA    \cr
##                           CA^2  \cr
##                           \vdots  \cr
##                           CA^{n-1} } \right ] $$
## @end tex
## @end iftex
## @ifinfo
## @example
## @group
##      | C        |
##      | CA       |
## Ob = | CA^2     |
##      | ...      |
##      | CA^(n-1) |
## @end group
## @end example
## @end ifinfo
## of a system data structure or the pair (@var{a}, @var{c}).
##
## The numerical properties of @command{is_observable}
## are much better for observability tests.
## @end deftypefn

## Author: Kai P. Mueller <mueller@ifr.ing.tu-bs.de>
## Created: November 4, 1997

function ob = obsv (sys_or_a, c)

  if (nargin == 1 && isstruct (sys))
    sysupdate (sys, "ss");
    [a, b, c] = sys2ss (sys);
  elseif (nargin == 2)
    a = sys_or_a;
    if (! isnumeric (a) || ! isnumeric (c) ||
        columns(a) != columns (c) || ! issquare (a))
      error ("obsv: invalid arguments");
    endif
  else
    print_usage ();
  endif

  ob = ctrb (a', c')';

endfunction


%!assert (obsv ([1 0; 0 -0.5], [8 8]), [8 8; 8 -4]);
