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
## @deftypefn{Function File} {@var{retsys} =} sysfeedback (@var{sys})
## Returns the closed loop of a system. The feedback is negative.
## Works for MIMO systems as well as for systems with mixed continuous
## and discrete parts. Replacement for unitfeedback.
##
## @strong{Inputs}
## @table @var
## @item sys
## System data structure. The number of inputs and outputs must be equal.
## @end table
##
## @strong{Outputs}
## @table @var
## @item retsys
## Closed loop system data structure.
## @end table
##
## @seealso{sysdup, sysscale, sysconnect, sysprune}
## 
## @example
## @group
##  u    +         +-------+         y
## ------>(+)----->|  sys  |-----+----->
##         ^ -     +-------+     |                             
##         |                     |
##         +---------------------+
##
##                      Y(s)       SYS(s)
##         RETSYS(s) = ------ = ------------
##                      U(s)     I + SYS(s)
##
##
## A SIMPLE EXAMPLE
##
##         P(s) : Plant          L(s) : Open Loop
##         C(s) : Controller     T(s) : Closed Loop
##
##  r    +    e    +--------+    u  +--------+         y
## ------>(+)----->|  C(s)  |------>|  P(s)  |-----+----->
##         ^ -     +--------+       +--------+     |
##         |                                       |
##         +---------------------------------------+
##
##         L(s) = P(s) C(s)      ---> L = sysmult(P, C)
##
##                   L(s)
##         T(s) = ----------     ---> T = sysfeedback(L)
##                 1 + L(s)
##
##  r      +--------+       y
## ------->|  T(s)  |--------->
##         +--------+
## @end group
## @end example
##
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@swissonline.ch>
## Version: 0.1


function retsys = sysfeedback (sys)

  if (nargin != 1)
    print_usage ();
  endif

  if(! isstruct (sys))
    error ("sysfeedback: argument sys must be a system data structure");
  endif

  [n_c_states, n_d_states, n_in, n_out] = sysdimensions (sys);
  
  if (n_in != n_out)
    error ("sysfeedback: number of inputs and outputs must be equal");
  endif
  
  ## Duplicate inputs
  for k = 1 : n_in
    in_idx(k) = k;
  endfor
  
  sys = sysdup (sys, [], in_idx);
  
  ## Scale inputs
  for k = 1 : n_in
    in_scl(k) = 1;
  endfor
  
  for k = (n_in + 1) : (2 * n_in)
    in_scl(k) = -1;
  endfor
  
  sys = sysscale (sys, [], diag (in_scl));
  
  ## Connect outputs with inputs
  for k = 1 : n_out
    out_idx(k) = k;
  endfor
  
  for k = 1 : n_in
    in_dup_idx(k) = k + n_in;
  endfor
  
  sys = sysconnect (sys, out_idx, in_dup_idx);
  
  ## Extract closed loop system
  retsys = sysprune (sys, out_idx, in_idx);
  
endfunction
