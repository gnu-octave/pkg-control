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
## @deftypefn {Function File} {@var{sys} =} sminreal (@var{sys})
## Perform state-space model reduction based on structure,
## i.e. retain only states which are both controllable and observable.
## The physical meaning of the states is retained.
## @end deftypefn

## Algorithm based on sysmin by A. Scottedward Hodel
## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.1

function sys = sminreal (sys)

  if (! isa (sys, "ss"))
    warning ("sminreal: system not in state-space form");
    sys = ss (sys);
  endif

  [cflg, Uc] = isctrb (sys);
  [oflg, Uo] = isobsv (sys);

  xc = find (max (abs (Uc.')) != 0);
  xo = find (max (abs (Uo.')) != 0);
  st_idx = intersect (xc, xo);

  sys = __sysprune__ (sys, ":", ":", st_idx);

  warning ("sminreal: use result with caution");

endfunction


## FIXME: algorithm returns wrong result for the example below
##
## P = ss (-2, 3, 4, 5)
## C = inv (P)
## L = P * C
## Ls = sminreal (L)
## bode (Ls)
## 
## Obviously, both magnitude and phase should be straight lines at 0 dB and 0 deg
## In this case, sminreal shouldn't remove any states of L