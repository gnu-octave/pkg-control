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
## Frequency response of FRD models :-)

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2010
## Version: 0.3

function H = __freqresp__ (sys, w, cellflag = false)

  [H, w_sys, tsam] = frdata (sys, "array");

  if (! isempty (w))     # freqresp (frdsys, w), sigma (frdsys, w), ...
    tol = sqrt (eps);
    w_idx = arrayfun (@(x) find (abs (w_sys - x) < tol), w, "uniformoutput", false);
    w_idx = vertcat (w_idx{:});

    ## NOTE: There are problems when cellfun uses "uniformoutput", true
    ##       and find returns an empty matrix,    

    if (length (w_idx) != numel (w))
      error ("frd: freqresp: some frequencies are not within tolerance %g", tol);
    endif

    H = H(:, :, w_idx);
  endif

  if (cellflag)
    [p, m, l] = size (H);
    H = mat2cell (H, p, m, ones (1, l))(:);
  endif

endfunction
