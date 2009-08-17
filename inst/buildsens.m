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
## @deftypefn{Function File} {@var{retsys} =} buildsens (@var{sys})
## Return the sensitivity of a system. Useful in combination with svplot.
##
## @strong{Inputs}
## @table @var
## @item sys
## Open loop system data structure. The number of inputs and outputs must be equal.
## @end table
##
## @strong{Outputs}
## @table @var
## @item retsys
## Sensitivity system data structure.
## @end table
##
## @seealso{svplot}
##
## @example
## @group
##         Y(s)         1               1
## S(s) = ------ = ---------- = -----------------
##         D(s)     1 + L(s)     1  +  P(s) C(s)
##
##
##
##                                                   | d
##                                                   |
##  r    +    e    +--------+    u  +--------+       v        y
## ------>(+)----->|  C(s)  |------>|  P(s)  |----->(+)---+----->
##         ^ -     +--------+       +--------+            |
##         |                                              |
##         +----------------------------------------------+
##
## L(s) = P(s) C(s)      ---> L = sysmult(P, C)
##
##            1
## S(s) = ----------     ---> S = buildsens(L)
##         1 + L(s)
##
##  d      +--------+       y
## ------->|  S(s)  |--------->
##         +--------+
## @end group
## @end example
##
## @end deftypefn

## Version: 0.1.1

function retsys = buildsens (sys)

  if (nargin != 1)
    print_usage ();
  endif

  if(! isstruct (sys))
    error ("buildsens: argument sys must be a system data structure");
  endif

  [n_c_states, n_d_states, n_in, n_out] = sysdimensions (sys);
  t_sam = sysgettsam (sys);
  digital = is_digital (sys, 2);

  if (n_in != n_out)
    error ("buildsens: number of inputs and outputs must be equal");
  endif

  if (digital == -1)
    error ("buildsens: system must be either purely continuous or purely discrete");
  endif

  sysI = ss ([], [], [], eye (n_in), t_sam);

  for k = 1 : n_in
    in_scl(k) = -1;
  endfor

  sys = sysscale (sys, [], diag (in_scl));

  retsys = feedback (sysI, sys, +1);

endfunction


%!shared numS, denS, numS_exp, denS_exp
%!
%! numL = [1];
%! denL = [1, 1];
%! L = tf(numL, denL);
%!
%! S = buildsens(L);
%! [numS, denS] = sys2tf(S);
%!
%! numS_exp = [1, 1];
%! denS_exp = [1, 2];
%!
%!assert(numS, numS_exp, eps);
%!assert(denS, denS_exp);
