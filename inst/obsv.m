## Copyright (C) 1997, 2000, 2002, 2004, 2005, 2006, 2007 Kai P. Mueller
## Copyright (C) 2009   Lukas F. Reichlin
## Copyright (C) 2009 Luca Favatella <slackydeb@gmail.com>
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
## @deftypefn {Function File} {@var{ob} =} obsv (@var{sys})
## @deftypefnx {Function File} {@var{ob} =} obsv (@var{a}, @var{c})
## Observability matrix.
## @end deftypefn

## Author: Kai P. Mueller <mueller@ifr.ing.tu-bs.de>
## Created: November 4, 1997

function ob = obsv (a, c)

  if (nargin == 1)  # obsv (sys)
    if (! isa (a, "lti"))
      error ("obsv: argument must be an lti system");
    endif
    [a, b, c] = ssdata (a);
  elseif (nargin == 2)  # obsv (a, c)
    if (! isnumeric (a) || ! isnumeric (c) ||
        columns (a) != columns (c) || ! issquare (a))
      error ("obsv: invalid arguments (a, c)");
    endif
  else
    print_usage ();
  endif

  ob = ctrb (a.', c.').';

endfunction


%!assert (obsv ([1, 0; 0, -0.5], [8, 8]), [8, 8; 8, -4]);
