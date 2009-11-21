## Copyright (C) 1997, 2000, 2002, 2004, 2005, 2006, 2007 Kai P. Mueller
## Copyright (C) 2009   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{co} =} ctrb (@var{sys})
## @deftypefnx {Function File} {@var{co} =} ctrb (@var{a}, @var{b})
## Build controllability matrix:
## @iftex
## @tex
## $$ C_o = [ B AB A^2B \ldots A^{n-1}B ] $$
## @end tex
## @end iftex
## @ifinfo
## @example
##              2       n-1
## Co = [ B AB A B ... A   B ]
## @end example
## @end ifinfo
##
## of a system data structure or the pair (@var{a}, @var{b}).
##
## @command{ctrb} forms the controllability matrix.
## The numerical properties of @command{is_controllable}
## are much better for controllability tests.
## @end deftypefn

## Author: Kai P. Mueller <mueller@ifr.ing.tu-bs.de>
## Created: November 4, 1997
## based on is_controllable.m of Scottedward Hodel

function co = ctrb (sys_or_a, b)

  if (nargin == 1 && isstruct (sys))
    sysupdate (sys, "ss");
    [a, b] = sys2ss (sys);
  elseif (nargin == 2)
    a = sys_or_a;
    if (! isnumeric (a) || ! isnumeric (b) ||
        rows(a) != rows (b) || ! issquare (a))
      error ("ctrb: invalid arguments");
    endif
  else
    print_usage ();
  endif

  [arows, acols] = size (a);
  [brows, bcols] = size (b);

  co = zeros (arows, acols*bcols);

  for k = 1 : arows
    co(:, ((k-1)*bcols + 1) : (k*bcols)) = b;
    b = a * b;
  endfor

endfunction


%!assert (ctrb ([1, 0; 0, -0.5], [8; 8]), [8, 8; 8, -4]);
