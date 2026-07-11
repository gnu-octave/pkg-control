## Copyright (C) 2026  Mitchell Thompkins <mitchell.thompkins@pm.me>
##
## This file is part of the control package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not,
## see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## Multiply two ZPK models element-wise.

## Author: Mitchell Thompkins <mitchell.thompkins@pm.me>
## Created: July 2026
## Version: 0.2

function sys = __times__ (sys1, sys2)
  if (! isa (sys1, "zpk"))
    sys1 = zpk (sys1);
  endif
  if (! isa (sys2, "zpk"))
    sys2 = zpk (sys2);
  endif

  [p1, m1] = size (sys1);
  [p2, m2] = size (sys2);

  if (p1 != p2 || m1 != m2)
    error ("zpk: __times__: system dimensions incompatible: (%dx%d) .* (%dx%d)",
            p1, m1, p2, m2);
  endif

  z = cellfun (@(a, b) [a; b], sys1.z, sys2.z, "uniformoutput", false);
  p = cellfun (@(a, b) [a; b], sys1.p, sys2.p, "uniformoutput", false);
  k = sys1.k .* sys2.k;

  ltisys = __lti_group__ (sys1.lti, sys2.lti, "times");

  sys = class (struct ("z", {z}, "p", {p}, "k", k), "zpk", ltisys);
endfunction
