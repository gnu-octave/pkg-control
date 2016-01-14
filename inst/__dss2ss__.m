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
## Convert descriptor state-space system into regular state-space form.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2011
## Version: 0.3

function [a, b, c, d, e] = __dss2ss__ (a, b, c, d, e)

  if (isempty (e))
    return;
  elseif (rcond (e) < eps)  # check for singularity
    ## check whether regular state-space representation is possible
    [~, ~, ~, ~, ranke, rnka22] = __sl_tg01fd__ (a, e, b, c, false, 0);
    if (ranke+rnka22 < rows (a))
      error ("dss:improper", "ss: dss2ss: this descriptor system cannot be converted to regular state-space form");
    endif
  endif

  [a, b, c, d] = __sl_sb10jd__ (a, b, c, d, e);
  e = [];

endfunction


## Test from SLICOT TG01FD
%!shared a, b, c, e, ranke, rnka22, q, z, a_exp, b_exp, c_exp, e_exp, q_exp, z_exp
%! 
%! e = [1, 2, 0, 0; 0, 1, 0, 1; 3, 9, 6, 3; 0, 0, 2, 0];
%! a = [-1, 0, 0, 3; 0, 0, 1, 2; 1, 1, 0, 4; 0, 0, 0, 0];
%! b = [1, 0; 0, 0; 0, 1; 1, 1];
%! c = [-1, 0, 1, 0; 0, 1, -1, 1];
%! 
%! [a, e, b, c, ranke, rnka22, q, z] = __sl_tg01fd__ (a, e, b, c, true, 0.0);
%! 
%! e_exp = [10.1587   5.8230   1.3021   0.0000;
%!           0.0000  -2.4684  -0.1896   0.0000;
%!           0.0000   0.0000   1.0338   0.0000;
%!           0.0000   0.0000   0.0000   0.0000];
%! 
%! a_exp = [ 2.0278   0.1078   3.9062  -2.1571;
%!          -0.0980   0.2544   1.6053  -0.1269;
%!           0.2713   0.7760  -0.3692  -0.4853;
%!           0.0690  -0.5669  -2.1974   0.3086];
%! 
%! b_exp = [-0.2157  -0.9705;
%!           0.3015   0.9516;
%!           0.7595   0.0991;
%!           1.1339   0.3780];
%! 
%! c_exp = [ 0.3651  -1.0000  -0.4472  -0.8165;
%!          -1.0954   1.0000  -0.8944   0.0000];
%! 
%! q_exp = [-0.2157  -0.5088   0.6109   0.5669;
%!          -0.1078  -0.2544  -0.7760   0.5669;
%!          -0.9705   0.1413  -0.0495  -0.1890;
%!           0.0000   0.8102   0.1486   0.5669];
%! 
%! z_exp = [-0.3651   0.0000   0.4472   0.8165;
%!          -0.9129   0.0000   0.0000  -0.4082;
%!           0.0000  -1.0000   0.0000   0.0000;
%!          -0.1826   0.0000  -0.8944   0.4082];
%! 
%!assert (a, a_exp, 1e-4);
%!assert (e, e_exp, 1e-4);
%!assert (b, b_exp, 1e-4);
%!assert (c, c_exp, 1e-4);
%!assert (q, q_exp, 1e-4);
%!assert (z, z_exp, 1e-4);
%!assert (ranke, 3);
%!assert (rnka22, 1);

## test error
%!shared mms
%! 
%! mm = tf ([3, 5, 0], [4, 1]);
%! mms = ss (mm);
%!error (__dss2ss__ (mms.a, mms.b, mms.c, mms.d, mms.e));

## Realizable descriptor system with singular E matrix
%!test
%! A = [1 0; 0 1];
%! B = [1; 0];
%! C = [1 0];
%! D = 0;
%! E = [1 0; 0 0];
%! 
%! sys = dss (A, B, C, D, E);
%! [Ao, Bo, Co, Do] = ssdata (sys);
%!
%! assert (Ao, 1, 1e-4);
%! assert (Bo, 1, 1e-4);
%! assert (Co, 1, 1e-4);
%! assert (Do, 0, 1e-4);
