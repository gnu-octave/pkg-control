## Copyright (C) 2021 Torsten Lilge <ttl-octave@mailbox.org>
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
## Make all corresponding num and den polynomials of a transfer function
## equally long by adding leading zeros to the shorter one.
## For internal use only.


function [num, den] = __make_tf_polys_equally_long__ (sys)

  num = sys.num;
  den = sys.den;

  for i = 1:size (num,1)
    for j = 1:size (num,2)
      [num{i,j},den{i,j}] = __make_equally_long__ (num{i,j},den{i,j});
    endfor
  endfor

endfunction
