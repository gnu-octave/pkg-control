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
## Block diagonal concatenation of two SS models.
## This file is part of the Model Abstraction Layer.
## For internal use only.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.1

function retsys = __sysgroup__ (sys1, sys2)

  if (! isa (sys1, "ss"))
    sys1 = ss (sys1);
  endif

  if (! isa (sys2, "ss"))
    sys2 = ss (sys2);
  endif

  retsys = ss ();

  retsys.lti = __ltigroup__ (sys1.lti, sys2.lti);

  A1 = sys1.a;
  B1 = sys1.b;
  C1 = sys1.c;
  D1 = sys1.d;
  A2 = sys2.a;
  B2 = sys2.b;
  C2 = sys2.c;
  D2 = sys2.d;

  [m1, n1, p1] = __ss_dim__ (A1, B1, C1, D1);
  [m2, n2, p2] = __ss_dim__ (A2, B2, C2, D2);

  A12 = zeros (n1, n2);
  B12 = zeros (n1, m2);
  B21 = zeros (n2, m1);
  C12 = zeros (p1, n2);
  C21 = zeros (p2, n1);
  D12 = zeros (p1, m2);
  D21 = zeros (p2, m1);

  retsys.a = [A1   , A12;
              A12.', A2 ];

  retsys.b = [B1 , B12;
              B21, B2 ];

  retsys.c = [C1 , C12;
              C21, C2 ];

  retsys.d = [D1 , D12;
              D21, D2 ];

  retsys.stname = [sys1.stname;
                   sys2.stname];

endfunction