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
## Frequency response of SS models.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.6

function H = __freqresp__ (sys, w, cellflag = false)

  if (sys.scaled == false)
    sys = prescale (sys);
  endif

  [a, b, c, d, e, tsam] = dssdata (sys);

  if (isct (sys))  # continuous system
    s = i * w;
  else             # discrete system
    s = exp (i * w * abs (tsam));
  endif

  ## compute H and prevent warnings when det(e*jw-A) gets close
  ## to a singular matrix near to resonance frequencies
  warning_id_near_sing = "Octave:nearly-singular-matrix";
  warning_near_sing = warning ("query", warning_id_near_sing);
  warning ("off", warning_id_near_sing);

  H = arrayfun (@(x) c/(x*e - a)*b + d, s, "uniformoutput", false);

  warning (warning_near_sing.state, warning_id_near_sing);

  if (! cellflag)
    H = cat (3, H{:});
  endif

endfunction
