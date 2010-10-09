## Copyright (C) 2010   Lukas F. Reichlin
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
## Block diagonal concatenation of two FRD models.
## This file is part of the Model Abstraction Layer.
## For internal use only.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2010
## Version: 0.1

function retsys = __sys_group__ (sys1, sys2)

  if (! isa (sys1, "frd"))
    sys1 = frd (sys1, sys2.w);
  endif

  if (! isa (sys2, "frd"))
    sys2 = frd (sys2, sys1.w);
  endif

  retsys = frd ();
  retsys.lti = __lti_group__ (sys1.lti, sys2.lti);

  lw1 = length (sys1.w);
  lw2 = length (sys2.w);

  [p1, m1, l1] = size (sys1.H);
  [p2, m2, l2] = size (sys2.H);


  if (lw1 == lw2 && all (sys1.w == sys2.w))
    retsys.w = sys1.w;
  else
    retsys.w = intersect (sys1.w, sys2.w);
  endif

  % l filtern
  H1 = sys1.H(:, :, w1_idx);
  H2 = sys2.H(:, :, w2_idx);

  z12 = zeros (p1, m2);
  z21 = zeros (p2, m1);

%  H1 = mat2cell (H1, p1, m1, ones (1, l))(:);
%  H2 = mat2cell (H2, p2, m2, ones (1, l))(:);
  H1 = mat2cell (H1, p1, m1, ones (1, l1))(:);
  H2 = mat2cell (H2, p2, m2, ones (1, l2))(:);

  H = cellfun (@(x, y) [x, z12; z21, y], H1, H2, "uniformoutput", false);

  retsys.H = cat (3, H{:});

% welche indices von H1, H2 ? --> vorher ausmisten

endfunction