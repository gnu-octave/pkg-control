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
## Block diagonal concatenation of two ZPK models.

## Author: Mitchell Thompkins <mitchell.thompkins@pm.me>
## Created: July 2026
## Version: 0.2

function retsys = __sys_group__ (sys1, sys2)
  [sys1, sys2] = __numeric_to_lti__ (sys1, sys2);

  if (! isa (sys1, "zpk"))
    sys1 = zpk (sys1);
  endif
  if (! isa (sys2, "zpk"))
    sys2 = zpk (sys2);
  endif

  [p1, m1] = size (sys1);
  [p2, m2] = size (sys2);

  z = [sys1.z, cell(p1, m2); cell(p2, m1), sys2.z];
  p = [sys1.p, cell(p1, m2); cell(p2, m1), sys2.p];
  k = [sys1.k, zeros(p1, m2); zeros(p2, m1), sys2.k];

  ltisys = __lti_group__ (sys1.lti, sys2.lti);

  retsys = class (struct ("z", {z}, "p", {p}, "k", k), "zpk", ltisys);
endfunction
