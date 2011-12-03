## Copyright (C) 2009, 2010, 2011   Lukas F. Reichlin
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
## Version: 0.3

function H = __freqresp__ (sys, w, resptype = 0, cellflag = false)

  [num, den, tsam] = tfdata (sys, "vector");

  if (isct (sys))  # continuous system
    s = num2cell (i * w);
  else             # discrete system
    s = num2cell (exp (i * w * abs (tsam)));
  endif

  if (issiso (sys))
    H = cellfun (@(z) polyval (num, z) / polyval (den, z), s, "uniformoutput", false);
  else
    f = @(z) cellfun (@(x, y) polyval (x, z) / polyval (y, z), num, den);
    H = cellfun (f, s, "uniformoutput", false);
  endif

  if (resptype)
    [p, m] = size (sys);

    if (m != p)
      error ("tf: freqresp: system must be square for response type %d", resptype);
    endif

    j = eye (p);

    switch (resptype)
      case 1     # inversed system
        H = cellfun (@inv, H, "uniformoutput", false);

      case 2     # inversed sensitivity
        H = cellfun (@(x) j + x, H, "uniformoutput", false);

      case 3     # inversed complementary sensitivity
        H = cellfun (@(x) j + inv (x), H, "uniformoutput", false);

      otherwise
        error ("tf: freqresp: invalid response type");
    endswitch
  endif

  if (! cellflag)
    H = cat (3, H{:});
  endif

endfunction
