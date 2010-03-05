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
## Frequency response of TF models.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.1

function H = __freqresp__ (sys, w, resptype = 0)

  [p, m] = size (sys);
  [num, den, Ts] = tfdata (sys);

  if (Ts > 0)  # discrete system
    s = exp (i * w * Ts);
  else  # continuous system
    s = i * w;
  endif

  l_s = length (s);
  H = zeros (p, m, l_s);

  for b = 1 : p
    for a = 1 : m
      num_pm = num{b, a};
      den_pm = den{b, a};

      for k = 1 : l_s
        H(b, a, k) = polyval (num_pm, s(k)) / polyval (den_pm, s(k));
      endfor
    endfor
  endfor

  if (resptype)
    if (m != p)
      error ("tf: freqresp: system must be square for response type %d", resptype);
    endif

    I = eye (p);

    switch (resptype)
      case 1  # inversed system
        for k = 1 : l_s
          H(:, :, k) = inv (H(:, :, k));
        endfor

      case 2  # inversed sensitivity
        for k = 1 : l_s
          H(:, :, k) = I + H(:, :, k);
        endfor

      case 3  # inversed complementary sensitivity
        for k = 1 : l_s
          H(:, :, k) = I + inv(H(:, :, k));
        endfor

      otherwise
        error ("tf: freqresp: invalid response type");

    endswitch
  endif

endfunction
