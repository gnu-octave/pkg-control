## Copyright (C) 2010   Lukas F. Reichlin
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
## Common code for adjusting FRD model data.
## Used by @frd/frd.m and @frd/__set__.m

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2010
## Version: 0.1

function [H, w, tsam] = __adjust_frd_data__ (H, w, tsam);

  w = reshape (w, [], 1);

  if (ndims (H) != 3 && ! isempty (H))
    if (is_real_scalar (H))           # static gain (H is a real scalar)
      H = reshape (H, 1, 1, []);
      tsam = -1;
    elseif (isvector (H))             # SISO system (H is a complex vector)
      H = reshape (H, 1, 1, []);
    elseif (is_real_matrix (H))       # static gain (H is a real matrix)
      H = reshape (H, rows (H), []);
      lw = length (w);
      if (lw > 1)
        H = repmat (H, [1, 1, lw]);   # needed for "frd1 + matrix2" or "matrix1 * frd2) 
      endif
      tsam = -1;
    else
      error ("frd: static gain matrix must be real");
    endif
  elseif (isempty (H))
    H = zeros (0, 0, 0);
    tsam = -1;
  endif

endfunction