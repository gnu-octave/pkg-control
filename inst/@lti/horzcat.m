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
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## Horizontal concatenation of LTI objects. If necessary, object conversion
## is done by sys_group. Used by Octave for "[lti1, lti2]".

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.1

function sys = horzcat (sys, varargin)

  for k = 1 : (nargin-1)

    sys1 = sys;
    sys2 = varargin{k};

    [p1, m1] = size (sys1);
    [p2, m2] = size (sys2);

    if (p1 != p2)
      error ("lti: horzcat: number of system outputs incompatible: [(%dx%d), (%dx%d)]",
              p1, m1, p2, m2);
    endif

    sys = __sys_group__ (sys1, sys2);

    out_scl = [eye(p1), eye(p2)];

    sys = out_scl * sys;

  endfor

endfunction