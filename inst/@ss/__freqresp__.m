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
## Frequency response of SS models.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.1

function H = __freqresp__ (sys, w, resptype = 0)

  [p, m] = size (sys);
  [A, B, C, D, Ts] = ssdata (sys);
  I = eye (size (A));
  J = eye (m);

  if (resptype != 0 && m != p)
    error ("ss: freqresp: system must be square for response type %d", resptype);
  endif

  if (Ts > 0)  # discrete system
    s = exp (i * w * Ts);
  else  # continuous system
    s = i * w;
  endif

  l_s = length (s);
  H = zeros (p, m, l_s);

  switch (resptype)
    case 0  # default system
      for k = 1 : l_s
        H(:, :, k) = C * inv (s(k)*I - A) * B  +  D;
      endfor

    case 1  # inversed system
      for k = 1 : l_s
        H(:, :, k) = inv (C * inv (s(k)*I - A) * B  +  D);
      endfor

    case 2  # inversed sensitivity
      for k = 1 : l_s
        H(:, :, k) = J  +  C * inv (s(k)*I - A) * B  +  D;
      endfor

    case 3  # inversed complementary sensitivity
      for k = 1 : l_s
        H(:, :, k) = J  +  inv (C * inv (s(k)*I - A) * B  +  D);
      endfor

    otherwise
      error ("ss: freqresp: invalid response type");

  endswitch

endfunction