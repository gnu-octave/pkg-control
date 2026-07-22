## Copyright (C) 2026        Prateek Ganguli
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

## Author: Prateek Ganguli <prateek.ganguli@gmail.com>
## Created: July 2026
## Version: 0.1

function y = __apply_timeresp_delay__ (y, delay)

  delay = round (delay);

  if (delay <= 0)
    return;
  endif

  n = rows (y);

  if (delay >= n)
    y = zeros (size (y));
  else
    y = [zeros(delay, columns (y)); y(1:end-delay, :)];
  endif

endfunction


%!test  # zero delay is a no-op
%! y = [1;2;3;4;5];
%! assert (__apply_timeresp_delay__ (y, 0), y);

%!test  # integer-sample shift, zero-padded at start, truncated at end
%! y = [1;2;3;4;5];
%! assert (__apply_timeresp_delay__ (y, 2), [0;0;1;2;3]);

%!test  # multi-column input shifts every column identically
%! y = [1 10; 2 20; 3 30; 4 40];
%! assert (__apply_timeresp_delay__ (y, 1), [0 0; 1 10; 2 20; 3 30]);

%!test  # delay >= number of rows yields all-zero output of the same size
%! y = [1;2;3];
%! assert (__apply_timeresp_delay__ (y, 5), zeros (3, 1));

%!test  # floating-point noise on an already-integer delay is tolerated
%! y = [1;2;3;4];
%! assert (__apply_timeresp_delay__ (y, 2 + 1e-10), [0;0;1;2]);
