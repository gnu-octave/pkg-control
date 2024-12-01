## Copyright (C) 2023  Torsten Lilge
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
## Chcking two systems that have to be connected if they
## possibly are only numerical values and have to be turned
## into proper lti systems.
## For internal use only.

## Author: Torsten Lilge <ttl-octave@mailbox.org>
## Created: July 2023
## Version: 0.1

function [sys1, sys2] = __numeric_to_lti__ (sys1, sys2)

  sys = {sys1, sys2};
  for i = 1:2
    if (! isa (sys{i}, "lti"))
      if (! isa (sys{i}, "numeric"))
        error ("lti: mtimes/mplus: one system is neither an lti system nor a numeric value");
      else
        sys{i} = tf (sys{i});
      endif
    endif
  endfor

  sys1 = sys{1};
  sys2 = sys{2};

  % If one of the two systems is only a static gain, just take the
  % sampling time of the other system
  if (isstaticgain(sys1) && (get(sys1,'tsam') == 0) && (get(sys2,'tsam') > 0))
    sys1 = c2d (sys1, get(sys2,'tsam'));
  elseif (isstaticgain(sys2)  && (get(sys2,'tsam') == 0) && (get(sys1,'tsam') > 0))
    sys2 = c2d (sys2, get(sys1,'tsam'));
  endif

endfunction

