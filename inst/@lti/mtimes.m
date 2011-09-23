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
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with LTI Syncope.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## Matrix multiplication of LTI objects.  If necessary, object conversion
## is done by sys_group.  Used by Octave for "sys1 * sys2".

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.1

function sys = mtimes (sys2, sys1)

  if (nargin != 2)  # prevent sys = mtimes (sys1, sys2, sys3, ...)
    error ("lti: mtimes: this is a binary operator");
  endif

  [p1, m1] = size (sys1);
  [p2, m2] = size (sys2);

  if (m2 != p1)
    error ("lti: mtimes: system dimensions incompatible: (%dx%d) * (%dx%d)",
            p2, m2, p1, m1);
  endif

  M22 = zeros (m2, p2);
  M21 = eye (m2, p1);
  M12 = zeros (m1, p2);
  M11 = zeros (m1, p1);

  M = [M22, M21;
       M12, M11];

  out_idx = 1 : p2;
  in_idx = m2 + (1 : m1);

  sys = __sys_group__ (sys2, sys1);
  sys = __sys_connect__ (sys, M);
  sys = __sys_prune__ (sys, out_idx, in_idx);

endfunction


## Alternative code: consistency vs. compatibility
#{
  M11 = zeros (m1, p1);
  M12 = zeros (m1, p2);
  M21 = eye (m2, p1);
  M22 = zeros (m2, p2);
  

  M = [M11, M12;
       M21, M22];

  out_idx = p1 + (1 : p2);
  in_idx = 1 : m1;

  sys = __sys_group__ (sys1, sys2);
#}
## Don't forget to adapt @tf/__sys_connect__.m draft code
