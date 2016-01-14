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
## Frequency response of TF models.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.5

function H = __freqresp__ (sys, w, cellflag = false)

  [num, den, tsam] = tfdata (sys, "vector");

  if (isct (sys))  # continuous system
    s = i * w;
  else             # discrete system
    s = exp (i * w * abs (tsam));
  endif
  
  s = reshape (s, 1, 1, []);

  if (issiso (sys))
    H = polyval (num, s) ./ polyval (den, s);
  else
    H = cellfun (@(x, y) polyval (x, s) ./ polyval (y, s), num, den, "uniformoutput", false);
    H = cell2mat (H);
  endif

  if (cellflag)
    [p, m] = size (sys);
    l = length (s);
    H = mat2cell (H, p, m, ones (1, l))(:);
  endif

endfunction
