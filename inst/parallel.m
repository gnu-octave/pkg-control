## Copyright (C) 2009 Lukas Reichlin <lukas.reichlin@swissonline.ch>
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn{Function File} {@var{sys} =} parallel (@var{sys1}, @var{sys2})
## Parallel connection of two systems.
## @math{sys = sys1 + sys2}
##
## @strong{Inputs}
## @table @var
## @item sys1
## System data structure or gain matrix.
## @item sys2
## System data structure or gain matrix.
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## System data structure of parallel connection.
## @end table
##
## @seealso{sysadd}
## 
## @example
## @group
##     ..........................
##     :      +--------+        :
##     :  +-->|  sys1  |---+    :
##  u  :  |   +--------+   | +  :  y
## -------+                O--------->
##     :  |   +--------+   | +  :
##     :  +-->|  sys2  |---+    :
##     :      +--------+        :
##     :.........sys............:
##
## sys = parallel (sys1, sys2)
## @end group
## @end example
##
## @end deftypefn

## Author: Lukas Reichlin
## Rewritten from scratch for better compatibility in July 2009
## Version: 0.1

function sys = parallel (_sys1, _sys2)

  if (nargin != 2)
    print_usage ();
  endif

  ## Sanitize sys1
  if (isstruct (_sys1))
    sys1 = _sys1;
    sys1wasmatrix = 0;
  elseif (ismatrix (_sys1))
    sys1 = ss ([], [], [], _sys1);
    sys1wasmatrix = 1;
  else
    error ("parallel: argument 1 (sys1) invalid");
  endif

  ## Sanitize sys2
  if (isstruct (_sys2))
    sys2 = _sys2;
    sys2wasmatrix = 0;
  elseif (ismatrix (_sys2))
    sys2 = ss ([], [], [], _sys2);
    sys2wasmatrix = 1;
  else
    error ("parallel: argument 2 (sys2) invalid");
  endif

  ## Handle digital gains
  if (sys1wasmatrix)
    if (is_digital (sys2, 2) == 1)  # -1 for mixed systems!
      t_sam = sysgettsam (sys2);
      sys1 = ss ([], [], [], _sys1, t_sam);  # sys1 = c2d (sys1, t_sam) doesn't work
    endif
  endif

  if (sys2wasmatrix)
    if (is_digital (sys1, 2) == 1)
      t_sam = sysgettsam (sys1);
      sys2 = ss ([], [], [], _sys2, t_sam);
    endif
  endif

  warning ("parallel: meaning has changed for compatibility reasons. use sysparallel for the old octave implementation.");

  ## Let sysadd do the job
  sys = sysadd (sys1, sys2);

endfunction


%!shared G, H, sys, num, den, num_exp, den_exp
%! G = ss ([0], [1], [1], [0]);
%! H = 1;
%! sys = parallel (G, H);
%! [num, den] = sys2tf (sys);
%! num_exp = [1, 1];
%! den_exp = [1, 0];
%!assert(num, num_exp);
%!assert(den, den_exp);
