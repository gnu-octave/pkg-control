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
## Inversion of SS models.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.1

function sys = __sys_inverse__ (sys)

  a = sys.a;
  b = sys.b;
  c = sys.c;
  d = sys.d;

  if (rcond (d) < eps)
    error ("ss: sys_inverse: inverse is not proper, case not implemented yet");
  else
    di = inv (d);

    f = a - b * di * c;
    g = b * di;
    h = -di * c;
    j = di;
  endif

  sys.a = f;
  sys.b = g;
  sys.c = h;
  sys.d = j;

endfunction