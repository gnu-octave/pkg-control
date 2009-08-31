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
## @deftypefn{Function File} {@var{sys} =} series (@var{sys1}, @var{sys2})
## @deftypefnx{Function File} {@var{sys} =} series (@var{sys1}, @var{sys2}, @var{outputs1}, @var{inputs2})
## Connect two systems in series.
## @math{sys = sys2 * sys1}
##
## @strong{Inputs}
## @table @var
## @item sys1
## System data structure or gain matrix.
## @item sys2
## System data structure or gain matrix.
## @item outputs1
## Optional vector of output indices for sys1.
## @item inputs2
## Optional vector of input indices for sys2.
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## System data structure of series connection.
## @end table
##
## @seealso{sysmult}
##
## @example
## @group
##     .....................................
##  u  :  +--------+ y1    u2  +--------+  :  y
## ------>|  sys1  |---------->|  sys2  |------->
##     :  +--------+           +--------+  :
##     :................sys.................
##
## sys = series (sys1, sys2)
##
##
##     .....................................
##     :                   v2  +--------+  :
##     :            ---------->|        |  :  y
##     :  +--------+ y1    u2  |  sys2  |------->
##  u  :  |        |---------->|        |  :
## ------>|  sys1  |       z1  +--------+  :
##     :  |        |---------->            :
##     :  +--------+                       :
##     :................sys.................
##
## outputs1 = [1]
## inputs2 = [2]
## sys = series (sys1, sys2, outputs1, inputs2)
## @end group
## @end example
##
## @end deftypefn

## Author: Lukas Reichlin
## Rewritten from scratch for better compatibility in July 2009
## Version: 0.2

function sys = series (_sys1, _sys2, outputs1, inputs2)

  if (nargin != 2 && nargin != 4)
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
    error ("series: argument 1 (sys1) invalid");
  endif

  ## Sanitize sys2
  if (isstruct (_sys2))
    sys2 = _sys2;
    sys2wasmatrix = 0;
  elseif (ismatrix (_sys2))
    sys2 = ss ([], [], [], _sys2);
    sys2wasmatrix = 1;
  else
    error ("series: argument 2 (sys2) invalid");
  endif

  ## Handle digital gains
  if (sys1wasmatrix)
    if (is_digital (sys2, 2) == 1) # -1 for mixed systems!
      t_sam = sysgettsam (sys2);
      sys1 = ss ([], [], [], _sys1, t_sam); # sys1 = c2d (sys1, t_sam) doesn't work
    endif
  endif

  if (sys2wasmatrix)
    if (is_digital (sys1, 2) == 1)
      t_sam = sysgettsam (sys1);
      sys2 = ss ([], [], [], _sys2, t_sam);
    endif
  endif


  if (nargin == 2)

    ## Let sysmult do the job
    sys = sysmult (sys2, sys1); # Note the different sequence of systems!

  else # (nargin == 4)

    ## Get system information
    [n_c_states_1, n_d_states_1, n_in_1, n_out_1] = sysdimensions (sys1);
    [n_c_states_2, n_d_states_2, n_in_2, n_out_2] = sysdimensions (sys2);

    ## Get connection lists outputs1 and inputs2
    if (! isvector (outputs1))
      error ("series: argument 3 (outputs1) invalid");
    endif

    if (! isvector (inputs2))
      error ("series: argument 4 (inputs2) invalid");
    endif

    l_outputs1 = length (outputs1);
    l_inputs2 = length (inputs2);

    if (l_outputs1 > n_out_1)
      error ("series: outputs1 has too many indices for sys1");
    endif

    if (l_inputs2 > n_in_2)
      error ("series: inputs2 has too many indices for sys2");
    endif

    if (l_outputs1 != l_inputs2)
      error ("series: number of outputs1 and inputs2 indices must be equal");
    endif


    ## Group sys1 and sys2 together

    ## Rename outputs of sys1 and inputs of sys2
    ## to avoid unnecessary warnings from sysgroup 

    out_name = __sysdefioname__ (n_out_1, "out1_series_tmp_name");
    in_name = __sysdefioname__ (n_in_2, "in2_series_tmp_name");

    sys1 = syssetsignals (sys1, "out", out_name);
    sys2 = syssetsignals (sys2, "in", in_name);

    sys = sysgroup (sys2, sys1); # Inversed sequence for compatibility


    ## Connect outputs of sys1 with inputs of sys2

    ## Output indices of sys1 start at (n_out_2 + 1)
    ## and end at (n_out_2 + n_out_1)
    ## Input indices of sys2 start at (1)
    ## and end at (n_in_2) --> use inputs2!

    out_idx = n_out_2 + outputs1;  

    sys = sysconnect (sys, out_idx, inputs2);


    ## Extract resulting model

    ## Output indices of sys2 start at (1)
    ## and end at (n_out_2)
    ## Input indices of sys1 start at (n_in_2 + 1)
    ## and end at (n_in_2 + n_in_1)

    out_pr_idx = 1 : n_out_2;

    in_pr_idx = (n_in_2 + 1) : (n_in_2 + n_in_1);

    sys = sysprune (sys, out_pr_idx, in_pr_idx);

  endif

endfunction


%!shared A1, B1, C1, D1, A2, B2, C2, D2, A_exp, B_exp, C_exp, D_exp
%! G = ss ([0, 1; -3, -2], [0; 1], [-5, 1], [2]);
%! H = ss ([-10], [1], [-40], [5]);
%! sys1 = series (G, H);
%! sys2 = series (G, H, 1, 1);
%! [A1, B1, C1, D1] = sys2ss (sys1);
%! [A2, B2, C2, D2] = sys2ss (sys2);
%! A_exp = [-10   -5    1 ;
%!            0    0    1 ;
%!            0   -3   -2 ];
%! B_exp = [  2 ;
%!            0 ;
%!            1 ];
%! C_exp = [-40  -25    5 ];
%! D_exp = [ 10 ];
%!assert([A1, B1; C1, D1], [A_exp, B_exp; C_exp, D_exp]);
%!assert([A2, B2; C2, D2], [A_exp, B_exp; C_exp, D_exp]);
