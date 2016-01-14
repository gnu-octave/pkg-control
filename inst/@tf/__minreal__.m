## Copyright (C) 2009-2016   Lukas F. Reichlin
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
## Minimal realization of TF models.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.2

function sys = __minreal__ (sys, tol)

  sqrt_eps = sqrt (eps);                        # treshold for zero
  [p, m] = size (sys);

  for ny = 1 : p
    for nu = 1 : m
      sisosys = __sys_prune__ (sys, ny, nu);

      [zer, gain] = zero (sisosys);
      pol = pole (sisosys);

      for k = length (zer) : -1 : 1             # reversed because of deleted zeros
        [~, idx] = min (abs (zer(k) - pol));    # find best match

        if (strcmpi (tol, "def"))
          if (abs (zer(k)) < sqrt_eps)          # catch case zer(k) = 0
            t = 1000 * eps;
          else
            t = 1000 * abs (zer(k)) * sqrt_eps;
          endif
        else
          t = tol;
        endif

        if (abs (zer(k) - pol(idx)) < t)
          zer(k) = [];
          pol(idx) = [];
        endif
      endfor

      num = real (gain * poly (zer));
      den = real (poly (pol));
      
      [n, d] = tfdata (sisosys, "vector");
      
      if (size_equal (num, n) && size_equal (den, d))
        if (d(1) != 1)
          n /= d(1);
          d /= d(1);
        endif
        num = n;
        den = d;
      endif

      sys.num{ny, nu} = tfpoly (num);
      sys.den{ny, nu} = tfpoly (den);
    endfor
  endfor

endfunction
