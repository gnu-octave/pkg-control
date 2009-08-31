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
## @deftypefn{Function File} {@var{retsys} =} buildretdiff (@var{sys})
## Build the return difference of a system.
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
## System data structure of the return difference.
## @end table
##
## @example
## @group
## Q(s) = I + L(s) = I  +  P(s) C(s)
## @end group
## @end example
##
## @end deftypefn

## Version: 0.1.3

function retsys = buildretdiff (sys)

  if (nargin != 1)
    print_usage ();
  endif

  if(! isstruct (sys))
    error ("buildretdiff: argument sys must be a system data structure");
  endif

  [n_c_states, n_d_states, n_in, n_out] = sysdimensions (sys);
  t_sam = sysgettsam (sys);
  digital = is_digital (sys, 2);

  if (n_in != n_out)
    error ("buildretdiff: number of inputs and outputs must be equal");
  endif

  if (digital == -1)
    error ("buildretdiff: system must be either purely continuous or purely discrete");
  endif

  sysI = ss ([], [], [], eye (n_in), t_sam);

  retsys = sysadd (sys, sysI);

endfunction


%!shared numQ, denQ, numQ_exp, denQ_exp
%!
%! numL = [1];
%! denL = [1, 1];
%! L = tf (numL, denL);
%!
%! Q = buildretdiff(L);
%! [numQ, denQ] = sys2tf (Q);
%!
%! numQ_exp = [1, 2];
%! denQ_exp = [1, 1];
%!
%!assert(numQ, numQ_exp, eps);
%!assert(denQ, denQ_exp);
