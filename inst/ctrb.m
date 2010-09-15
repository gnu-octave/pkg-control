## Copyright (C) 2009 - 2010   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{co} =} ctrb (@var{sys})
## @deftypefnx {Function File} {@var{co} =} ctrb (@var{a}, @var{b})
## Controllability matrix.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.2

function co = ctrb (a, b)

  if (nargin == 1)       # ctrb (sys)
    if (! isa (a, "lti"))
      error ("ctrb: argument must be an lti system");
    endif
    [a, b] = ssdata (a);
  elseif (nargin == 2)   # ctrb (a, b)
    if (! isreal (a) || ! isreal (b)
        || rows (a) != rows (b) || ! issquare (a))
      error ("ctrb: invalid arguments (a, b)");
    endif
  else
    print_usage ();
  endif

  n = rows (a);          # number of states
  k = num2cell (0:n-1);  # exponents for a

  tmp = cellfun (@(x) a^x*b, k, "uniformoutput", false);

  co = horzcat (tmp{:});

endfunction


%!assert (ctrb ([1, 0; 0, -0.5], [8; 8]), [8, 8; 8, -4]);
