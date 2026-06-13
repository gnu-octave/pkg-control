## Copyright (C) 2026 Dmitri A. Sergatskov <dasergatskov@gmail.com>
##
## This function is part of the GNU Octave Control Package
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## DC gain of TF models.

function gain = dcgain (sys)

  if (nargin != 1)
    print_usage ();
  endif

  [num, den] = tfdata (sys);

  if (isct (sys))
    gain = cellfun (@(n, d) n(end) ./ d(end), num, den);
  else
    gain = cellfun (@(n, d) polyval (n, 1) ./ polyval (d, 1), num, den);
  endif

endfunction

%!assert (dcgain (tf (1, [1, 1])), 1)
%!assert (dcgain (tf (2, [1, 1])), 2)
%!assert (dcgain (tf (1, [1, -0.5], 1)), 2)
%!assert (dcgain (tf (1, [1, -1], 1)), Inf)

%!test
%! num = [0, 0, 2.052877715426585e9, 4.416784109977014e13, 1.719831064785571e17];
%! den = [1, 8.380906856780281e4, 2.269817212624148e8, 4.344755797517798e12, 0];
%! gain = dcgain (tf (num, den));
%! assert (isinf (gain));
%! assert (gain > 0);
