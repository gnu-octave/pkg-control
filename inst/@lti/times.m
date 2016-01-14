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
## Hadamard/Schur product of transfer function matrices.
## Also known as element-wise multiplication.
## Used by Octave for "sys1 .* sys2".@*
## @strong{Example}
## @example
## @group
## # Compute Relative-Gain Array
## G = tf (Boeing707)
## RGA = G .* inv (G).'
## # Gain at 0 rad/s
## RGA(0)
## @end group
## @end example

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: April 2014
## Version: 0.2

function sys = times (sys1, sys2)

  if (nargin != 2)                              # prevent sys = times (sys1, sys2, sys3, ...)
    error ("lti: times: this is a binary operator");
  endif

  [p1, m1] = size (sys1);
  [p2, m2] = size (sys2);
  
  if (p1 != p2 || m1 != m2)
    if (p1 == 1 && m1 == 1 && p2*m2 > 1)        # sys1 SISO, sys2 non-empty
      sys1 = repmat (sys1, p2, m2);
    elseif (p2 == 1 && m2 == 1 && p1*m1 > 1)    # sys2 SISO, sys1 non-empty
      sys2 = repmat (sys2, p1, m1);
    else
      error ("lti: times: system dimensions incompatible: (%dx%d) .* (%dx%d)", ...
              p1, m1, p2, m2);
    endif
  endif
  
  sys = __times__ (sys1, sys2);

endfunction
