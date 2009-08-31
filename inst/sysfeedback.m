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
## and discrete parts.
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
##              Y(s)       SYS(s)
## RETSYS(s) = ------ = ------------
##              U(s)     I + SYS(s)
##
##
## A SIMPLE EXAMPLE
##
## P(s) : Plant          L(s) : Open Loop
## C(s) : Controller     T(s) : Closed Loop
##
##  r    +    e    +--------+    u  +--------+         y
## ------>(+)----->|  C(s)  |------>|  P(s)  |-----+----->
##         ^ -     +--------+       +--------+     |
##         |                                       |
##         +---------------------------------------+
##
## L(s) = P(s) C(s)      ---> L = sysmult(P, C)
##
##           L(s)
## T(s) = ----------     ---> T = sysfeedback(L)
##         1 + L(s)
##
##  r      +--------+       y
## ------->|  T(s)  |--------->
##         +--------+
## @end group
## @end example
##
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@swissonline.ch>
## Version: 0.2.1

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
  in_idx = 1 : n_in;

  sys = sysdup (sys, [], in_idx);

  ## Scale inputs
  in_scl = ones (1, n_in);
  in_scl = horzcat (in_scl, (-1) * in_scl);

  sys = sysscale (sys, [], diag (in_scl));

  ## Connect outputs with inputs
  out_idx = 1 : n_out;

  in_dup_idx = n_in + (1 : n_in);

  sys = sysconnect (sys, out_idx, in_dup_idx);

  ## Extract closed loop system
  retsys = sysprune (sys, out_idx, in_idx);

endfunction


%!test
%!
%! ## Given Open Loop State Space Matrices
%!
%! A = [        -25         0         0         0        50         0       763     310.7    -27.93    -96.51     216.4       -61 ;
%!            0.261        -4    -0.002         0         0         0         0         0         0         0         0         0 ;
%!              618      7640     -3.07         0         0     338.4     17.19     -2247      5459     -2782     298.7      2849 ;
%!           -0.384      5.89     0.003      -2.1         0     -0.21  -0.01822    0.9135   -0.9346     2.552   -0.0403    -1.761 ;
%!                0         0         0         0         0         0     30.52     12.43    -1.117    -3.861     8.658     -2.44 ;
%!                0         0         0         0         0         0    0.1735      -8.7     8.901    -24.31    0.3839     16.77 ;
%!                0         0         0         0         0         0      -788    -310.7     26.59     96.98    -166.4        61 ;
%!                0         0         0         0         0         0      5.22        -4    -7.945    0.2199         0         0 ;
%!                0         0         0         0         0         0     3.004     13.14    -31.54     13.91    -1.493    -12.56 ;
%!                0         0         0         0         0         0    -7.316    -12.38     30.69    -54.82    0.8061     31.01 ;
%!                0         0         0         0         0         0    -30.52    -12.43    0.6207     3.803    -8.658      2.44 ;
%!                0         0         0         0         0         0   -0.1735       8.7    -8.958      24.8   -0.3839    -16.77 ];
%!
%! B = [            0           0 ;
%!                  0           0 ;
%!                  0           0 ;
%!                  0           0 ;
%!                  0           0 ;
%!                  0           0 ;
%!           0.006741      -9.319 ;
%!         -0.0002757      -4.397 ;
%!           0.005862     0.01201 ;
%!          3.002e-06       33.54 ;
%!           0.002484       1.146 ;
%!          0.0002864      -9.934 ];
%!
%! C = [  0    0    1    0    0    0    0    0    0    0    0    0 ;
%!        0    0    0    1    0    0    0    0    0    0    0    0 ];
%!
%! D = [  0   0 ;
%!        0   0 ];
%!
%! ## Expected Closed Loop State Space Matrices 
%!
%! F_exp = [      -25           0           0           0          50           0         763       310.7      -27.93      -96.51       216.4         -61 ;
%!              0.261          -4      -0.002           0           0           0           0           0           0           0           0           0 ;
%!                618        7640       -3.07           0           0       338.4       17.19       -2247        5459       -2782       298.7        2849 ;
%!             -0.384        5.89       0.003        -2.1           0       -0.21    -0.01822      0.9135     -0.9346       2.552     -0.0403      -1.761 ;
%!                  0           0           0           0           0           0       30.52       12.43      -1.117      -3.861       8.658       -2.44 ;
%!                  0           0           0           0           0           0      0.1735        -8.7       8.901      -24.31      0.3839       16.77 ;
%!                  0           0   -0.006741       9.319           0           0        -788      -310.7       26.59       96.98      -166.4          61 ;
%!                  0           0   0.0002757       4.397           0           0        5.22          -4      -7.945      0.2199           0           0 ;
%!                  0           0   -0.005862    -0.01201           0           0       3.004       13.14      -31.54       13.91      -1.493      -12.56 ;
%!                  0           0  -3.002e-06      -33.54           0           0      -7.316      -12.38       30.69      -54.82      0.8061       31.01 ;
%!                  0           0   -0.002484      -1.146           0           0      -30.52      -12.43      0.6207       3.803      -8.658        2.44 ;
%!                  0           0  -0.0002864       9.934           0           0     -0.1735         8.7      -8.958        24.8     -0.3839      -16.77 ];
%!
%! G_exp = [        0           0 ;
%!                  0           0 ;
%!                  0           0 ;
%!                  0           0 ;
%!                  0           0 ;
%!                  0           0 ;
%!           0.006741      -9.319 ;
%!         -0.0002757      -4.397 ;
%!           0.005862     0.01201 ;
%!          3.002e-06       33.54 ;
%!           0.002484       1.146 ;
%!          0.0002864      -9.934 ];
%!
%! H_exp = [   0    0    1    0    0    0    0    0    0    0    0    0 ;
%!             0    0    0    1    0    0    0    0    0    0    0    0 ];
%!
%! J_exp = [   0   0 ;
%!             0   0 ];
%!
%! sys_ol = ss (A, B, C, D);
%! sys_cl = sysfeedback (sys_ol);
%! [F, G, H, J] = sys2ss (sys_cl);
%!
%! assert([F, G; H, J], [F_exp, G_exp; H_exp, J_exp]);
