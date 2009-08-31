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
## Return model sys for the negative feedback interconnection; i.e. filter the output of
## @var{sys1} through @var{sys2} and subtract it from the input.
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
## @code{-1} means negative feedback (default if not specified),
## @code{+1} means positive one.
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

## Author: Lukas Reichlin
## Rewritten from scratch for better compatibility in July 2009
## Version: 0.2.6

function sys = feedback (_sys1, _sys2, _sign_or_feedin, _feedout, sign = -1)

  if (nargin < 2 || nargin > 5)
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
    error ("feedback: argument 1 (sys1) invalid");
  endif

  ## Sanitize sys2
  if (isstruct (_sys2))
    sys2 = _sys2;
    sys2wasmatrix = 0;
  elseif (ismatrix (_sys2))
    sys2 = ss ([], [], [], _sys2);
    sys2wasmatrix = 1;
  else
    error ("feedback: argument 2 (sys2) invalid");
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

  ## Determine feedback sign
  if (nargin == 3)
    sign = _sign_or_feedin;
  endif

  ## Sanitize feedback sign
  if ((! isreal (sign)) && ((sign != -1) && # "&&" is short-circuit
                            (sign != +1)))
    error ("feedback: sign must be -1 or +1");
  endif

  ## Get system information
  [n_c_states_1, n_d_states, n_in_1, n_out_1] = sysdimensions (sys1);
  [n_c_states_2, n_d_states_2, n_in_2, n_out_2] = sysdimensions (sys2);


  ## Get connection lists feedin and feedout
  feedin = 1 : n_in_1; # Default value
  feedout = 1 : n_out_1; # Default value

  if (nargin == 4 || nargin == 5)
    if (isvector (_sign_or_feedin))
      feedin = _sign_or_feedin;
    else
      error ("feedback: argument 3 (feedin) invalid");
    endif

    if (isvector (_feedout))
      feedin = _feedout;
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

  ## Rename inputs, outputs and states of sys2 temporarily
  ## to avoid spurious warnings from sysgroup

  in_name = __sysdefioname__ (n_in_2, "in_feedback_tmp_name");
  out_name = __sysdefioname__ (n_out_2, "out_feedback_tmp_name");

  sys2 = syssetsignals (sys2, "in", in_name);
  sys2 = syssetsignals (sys2, "out", out_name);

  if ((n_c_states_2 + n_d_states_2) != 0) # If there are any states

    for k = 1 : (n_c_states_2 + n_d_states_2)
      st_name{k} = sprintf ("st_feedback_tmp_name_%d", k);
    endfor

    sys2 = syssetsignals (sys2, "st", st_name);

  endif

  sys = sysgroup (sys1, sys2);


  ## Duplicate inputs specified in feedin
  sys = sysdup(sys, [], feedin);

  ## Duplicated inputs start now at number (n_in_1 + n_in_2 + 1)
  ## and end at number (n_in_1 + n_in_2 + l_feedin)


  ## Scale inputs
  in_scl(1, 1 : (n_in_1 + n_in_2 + l_feedin)) = \
      [(ones (1, n_in_1 + n_in_2)), (sign * ones (1, l_feedin))];
  sys = sysscale (sys, [], diag (in_scl));


  ## Connect outputs with inputs
  ## FIXME: Unable to prevent algebraic loop warnings from sysconnect
  ## even if connections are made one by one with
  ## for k = 1 : l_feedout
  ##     sys = sysconnect (sys, feedout(k), in_idx(k));
  ## endfor
  ## and
  ## for k = 1 : l_feedin
  ##   sys = sysconnect (sys, out_idx(k), in_dup_idx(k));
  ## endfor

  ## Connect outputs of sys1 to inputs of sys2
  ## Output indices of sys1 start at (1)
  ## and end at (n_out_1) --> use feedout!
  ## Input indices of sys2 start at (n_in_1 + 1)
  ## and end at (n_in_1 + n_in_2)

  in_idx = n_in_1 + (1 : l_feedout);
  sys = sysconnect (sys, feedout, in_idx);

  ## Connect outputs of sys2 to inputs of sys1
  ## Output indices of sys 2 start at (n_out_1 + 1)
  ## and end at (n_out_1 + n_out_2)
  ## Duplicated input indices of sys 1 start at (n_in_1 + n_in_2 + 1)
  ## and end at (n_in_1 + n_in_2 + l_feedin)
  ## Sequence of feedin already mapped to duplicated inputs
  ## ---> (n_in_1 + n_in_2 + k)

  out_idx = n_out_1 + (1 : l_feedin);
  in_dup_idx = n_in_1 + n_in_2 + (1 : l_feedin);
  sys = sysconnect (sys, out_idx, in_dup_idx);


  ## Extract resulting model
  out_pr_idx = 1 : n_out_1;
  in_pr_idx = 1 : n_in_1;
  sys = sysprune (sys, out_pr_idx, in_pr_idx);

endfunction


%!shared sys1, sys2, sys3, sys4, sys5, sys6, sysexpa, sysexpb
%! G = ss ([0], [1], [1], [0]);
%! H = ss ([], [], [], [1]);
%!
%! sys1 = feedback (G, H);
%! sys2 = feedback (G, 1, -1);
%! sys3 = feedback (G, H, 1, 1, -1);
%! sysexpa = ss ([-1], [1], [1], [0]);
%!
%! sys4 = feedback (G, H, +1);
%! sys5 = feedback (G, 1, +1);
%! sys6 = feedback (G, H, 1, 1, +1);
%! sysexpb = ss ([1], [1], [1], [0]);
%!
%!assert (sys1, sysexpa);
%!assert (sys2, sysexpa);
%!assert (sys3, sysexpa);
%!assert (sys4, sysexpb);
%!assert (sys5, sysexpb);
%!assert (sys6, sysexpb);
