## Copyright (C) 2000 Ben Sapp <bsapp@lanl.gov>
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
## @deftypefn{Function File} {@var{sys} =} feedback (@var{sys1}, @var{sys2})
## @deftypefnx{Function File} {@var{sys} =} feedback (@var{sys1}, @var{sys2}, @var{sign})
## @deftypefnx{Function File} {@var{sys} =} feedback (@var{sys1}, @var{sys2}, @var{feedin}, @var{feedout})
## @deftypefnx{Function File} {@var{sys} =} feedback (@var{sys1}, @var{sys2}, @var{feedin}, @var{feedout}, @var{sign})
## Returns model sys for the feedback interconnection; i. e. filters the output of sys1
## through sys2 and subtracts it from the input. Displays warnings if inputs, outputs
## or states of sys1 and sys2 have equal names. Occasionally warns about possible algebraic
## loops as well. These warnings are perfectly harmless and are due to the internal use
## of sysgroup and sysconnect.
##
## @strong{Inputs}
## @table @var
## @item sys1
## System data structure or gain matrix.
## @item sys2
## System data structure or gain matrix.
## @item feedin
## Optional vector of input indices for sys1.
## @item feedout
## Optional vector of output indices for sys1.
## @item sign
## If not specified, default value -1 is taken, resulting in a negative feedback.
## Use +1 for positive feedback.
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## System data structure of feedback interconnection.
## @end table
##
## @seealso{sysgroup, sysdup, sysscale, sysconnect, sysprune}
## 
## @example
## @group
##  u    +         +--------+             y
## ------>(+)----->|  sys1  |-------+------->
##         ^ -     +--------+       |                             
##         |                        |
##         |       +--------+       |
##         +-------|  sys2  |<------+
##                 +--------+                             
## @end group
## @end example
##
## @end deftypefn
 
## Author: Ben Sapp <bsapp@lanl.gov>
## Rewritten from scratch by Lukas Reichlin
## for better M*tl*b compatibility in July 2009
## Version: 0.1
 
function sys = feedback (varargin)
  
  if (nargin < 2 || nargin > 5)
    print_usage ();
  endif
  
  ## Determine sys1 from input
  if (isstruct (varargin{1}))
    sys1 = varargin{1};
  elseif (ismatrix (varargin{1}))
    sys1 = ss ([], [], [], varargin{1});
  else
    error ("feedback: argument 1 (sys1) invalid");
  endif
  
  ## Determine sys2 from input
  if (isstruct (varargin{2}))
    sys2 = varargin{2};
  elseif (ismatrix (varargin{2}))
    sys2 = ss ([], [], [], varargin{2});
  else
    error ("feedback: argument 2 (sys2) invalid");
  endif
  
  ## Determine feedback sign
  fb_sign = -1; # Default value
  
  if (nargin == 3)
    if (isreal (varargin{3}))
      if (varargin{3} == +1)
        fb_sign = +1;
      endif
      if (varargin{3} != -1 && varargin{3} != +1)
        error ("feedback: argument 3 (sign) invalid");
      endif
    endif
  endif
  
  if (nargin == 5)
    if (isreal (varargin{5}))
      if (varargin{5} == +1)
        fb_sign = +1;
      endif
      if (varargin{5} != -1 && varargin{3} != +1)
        error ("feedback: argument 5 (sign) invalid");
      endif
    endif
  endif
  
  ## Get system information
  [n_c_states_1, n_d_states, n_in_1, n_out_1] = sysdimensions (sys1);
  [n_c_states_2, n_d_states_2, n_in_2, n_out_2] = sysdimensions (sys2);
  
  
  ## Get connection lists feedin and feedout
  for k = 1 : n_in_1
    feedin(k) = k; # Default value
  endfor
  
  for k = 1 : n_out_1
    feedout(k) = k; # Default value
  endfor
  
  if (nargin == 4 || nargin == 5)
    if (isvector (varargin{3}))
      feedin = varargin{3};
    else
      error ("feedback: argument 3 (feedin) invalid");
    endif
    
    if (isvector (varargin{4}))
      feedin = varargin{4};
    else
      error ("feedback: argument 4 (feedout) invalid");
    endif
  endif
  
  l_feedin = length (feedin);
  l_feedout = length (feedout);
  
  if (l_feedin > n_in_1)
    error ("feedback: feedin has too many indices");
  endif
  
  if (l_feedin != n_out_2)
    error ("feedback: number of feedin indices and number of outputs of sys2 incompatible");
  endif
  
  if (l_feedout > n_out_1)
    error ("feedback: feedout has too many indices");
  endif
  
  if (l_feedout != n_in_2)
    error ("feedback: number of feedout indices and number of inputs of sys2 incompatible");
  endif
  
  
  ## Group sys1 and sys2 together
  sys = sysgroup(sys1, sys2);
  
  
  ## Duplicate inputs specified in feedin
  sys = sysdup(sys, [], feedin);
  
  ## Duplicated inputs start now at number (n_in_1 + n_in_2 + 1)
  ## and end at number (n_in_1 + n_in_2 + l_feedin)
  
  
  ## Scale inputs
  for k = 1 : (n_in_1 + n_in_2)
    in_scl(k) = 1;
  endfor
  
  for k = (n_in_1 + n_in_2 + 1) : (n_in_1 + n_in_2 + l_feedin)
    in_scl(k) = fb_sign;
  endfor
   
  sys = sysscale (sys, [], diag (in_scl));
  
  
  ## Connect outputs with inputs
  
  ## connect outputs of sys1 to inputs of sys2
  ## output indices of sys1 start at (1)
  ## and end at (n_out_1) --> use feedout!
  ## input indices of sys2 start at (n_in_1 + 1)
  ## and end at (n_in_1 + n_in_2)
  
  for k = 1 : l_feedout
    in_idx(k) = n_in_1 + k;
  endfor
  
  sys = sysconnect (sys, feedout, in_idx);
  
  ## connect outputs of sys2 to inputs of sys1
  ## output indices of sys 2 start at (n_out_1 + 1)
  ## and end at (n_out_1 + n_out_2)
  ## duplicated input indices of sys 1 start at (n_in_1 + n_in_2 + 1)
  ## and end at (n_in_1 + n_in_2 + l_feedin)
  ## sequence of feedin already mapped to duplicated inputs
  ## ---> (n_in_1 + n_in_2 + k)

  for k = 1 : l_feedin
    out_idx(k) = n_out_1 + k;
  endfor
  
  for k = 1 : l_feedin
    in_dup_idx(k) = n_in_1 + n_in_2 + k;
  endfor
  
  sys = sysconnect (sys, out_idx, in_dup_idx);
  
  
  ## Extract resulting model
  for k = 1 : n_out_1
    out_pr_idx(k) = k;
  endfor
  
  for k = 1 : n_in_1
    in_pr_idx(k) = k;
  endfor
  
  sys = sysprune(sys, out_pr_idx, in_pr_idx);
  
endfunction


%!shared A, B, C, D, F_exp, G_exp, H_exp, J_exp, F, G, H, J, sys, sysc, M, N, outsys, zer, pol, k, zer_exp, pol_exp, k_exp
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
%! sys1 = ss(A, B, C, D);
%! sys = feedback(sys1, eye(2));
%! [F, G, H, J] = sys2ss(sys);
%! 
%! sysc = feedback(sys1, eye(2), [1, 2], [1, 2], -1);
%!
%!
%! M = tf([2, 5, 1], [1, 2, 3]);
%! N = zp(-2, -10, 5);
%! outsys = feedback(M, N);
%! [zer, pol, k] = sys2zp(outsys);
%!
%! zer_exp = [ -2.280776406404418 ;
%!             -0.219223593595584 ;
%!             -9.999999999999996 ];
%! pol_exp = [ -0.881474800310081 + 0.535367413587643i ;
%!             -0.881474800310081 - 0.535367413587643i ;
%!             -3.418868581198019 + 0.000000000000000i ];
%! k_exp =  0.181818181818182;
%!
%!assert(F, F_exp);
%!assert(G, G_exp);
%!assert(H, H_exp);
%!assert(J, J_exp);
%!assert(sys, sysc)
%!assert(zer, zer_exp);
%!assert(pol, pol_exp, 2*eps);
%!assert(k, k_exp, 2*eps);
