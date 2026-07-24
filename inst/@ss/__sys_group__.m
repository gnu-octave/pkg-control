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
## Block diagonal concatenation of two SS models.
## This file is part of the Model Abstraction Layer.
## For internal use only.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.2

function retsys = __sys_group__ (sys1, sys2)

  % If one system is just a numeric value, create a proper lti system
  [sys1, sys2] = __numeric_to_lti__ (sys1, sys2);

  if (! isa (sys1, "ss"))
    sys1 = ss (sys1);
  endif

  if (! isa (sys2, "ss"))
    sys2 = ss (sys2);
  endif

  retsys = ss ();
  retsys.lti = __lti_group__ (sys1.lti, sys2.lti);

  n1 = rows (sys1.a);
  n2 = rows (sys2.a);
  
  [p1, m1] = size (sys1.d);
  [p2, m2] = size (sys2.d);

  retsys.a = [sys1.a, zeros(n1,n2); zeros(n2,n1), sys2.a];
  retsys.b = [sys1.b, zeros(n1,m2); zeros(n2,m1), sys2.b];
  retsys.c = [sys1.c, zeros(p1,n2); zeros(p2,n1), sys2.c];
  retsys.d = [sys1.d, zeros(p1,m2); zeros(p2,m1), sys2.d];

  e1 = ! isempty (sys1.e);
  e2 = ! isempty (sys2.e);

  if (e1 || e2)
    if (e1 && e2)
      retsys.e = [sys1.e, zeros(n1,n2); zeros(n2,n1), sys2.e];
    elseif (e1)
      retsys.e = [sys1.e, zeros(n1,n2); zeros(n2,n1), eye(n2)];
    else
      retsys.e = [eye(n1), zeros(n1,n2); zeros(n2,n1), sys2.e];
    endif
  endif

  retsys.stname = [sys1.stname; sys2.stname];

  id1 = ! isempty (sys1.b2);
  id2 = ! isempty (sys2.b2);

  if (id1 || id2)
    if (id1)
      q1 = columns (sys1.b2);
      b21 = sys1.b2;  c21 = sys1.c2;
      d121 = sys1.d12;  d211 = sys1.d21;  d221 = sys1.d22;
    else
      q1 = 0;
      b21 = zeros (n1, 0);  c21 = zeros (0, n1);
      d121 = zeros (p1, 0);  d211 = zeros (0, m1);  d221 = zeros (0, 0);
    endif

    if (id2)
      q2 = columns (sys2.b2);
      b22 = sys2.b2;  c22 = sys2.c2;
      d122 = sys2.d12;  d212 = sys2.d21;  d222 = sys2.d22;
    else
      q2 = 0;
      b22 = zeros (n2, 0);  c22 = zeros (0, n2);
      d122 = zeros (p2, 0);  d212 = zeros (0, m2);  d222 = zeros (0, 0);
    endif

    ## delay ports are independent between operands: no cross-terms
    retsys.b2 = [b21, zeros(n1,q2); zeros(n2,q1), b22];
    retsys.c2 = [c21, zeros(q1,n2); zeros(q2,n1), c22];
    retsys.d12 = [d121, zeros(p1,q2); zeros(p2,q1), d122];
    retsys.d21 = [d211, zeros(q1,m2); zeros(q2,m1), d212];
    retsys.d22 = [d221, zeros(q1,q2); zeros(q2,q1), d222];
  endif

  retsys.internaldelay = [sys1.internaldelay; sys2.internaldelay];

endfunction
