## Copyright (C) 2009, 2011   Lukas F. Reichlin
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
## Convert the continuous TF model into its discrete-time equivalent.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.2

function sys = __c2d__ (sys, tsam, method = "zoh", w0 = 0)

  [p, m] = size (sys);

  for i = 1 : p
    for j = 1 : m
      idx = substruct ("()", {i, j});
      tmp = subsref (sys, idx);
      tmp = c2d (ss (tmp), tsam, method, w0);
      [num, den] = tfdata (tmp, "tfpoly");
      sys.num(i, j) = num;
      sys.den(i, j) = den;
    endfor
  endfor

  sys.tfvar = "z";

endfunction
