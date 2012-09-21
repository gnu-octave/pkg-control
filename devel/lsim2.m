## Copyright (C) 2009, 2010, 2011   Lukas F. Reichlin
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
## @deftypefn{Function File} {[@var{y}, @var{t}, @var{x}] =} lsim (@var{sys}, @var{u})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} lsim (@var{sys}, @var{u}, @var{t})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} lsim (@var{sys}, @var{u}, @var{t}, @var{x0})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} lsim (@var{sys}, @var{u}, @var{t}, @var{[]}, @var{method})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} lsim (@var{sys}, @var{u}, @var{t}, @var{x0}, @var{method})
## Simulate LTI model response to arbitrary inputs.  If no output arguments are given,
## the system response is plotted on the screen.
##
## @strong{Inputs}
## @table @var
## @item sys
## LTI model.  System must be proper, i.e. it must not have more zeros than poles.
## @item u
## Vector or array of input signal.  Needs @code{length(t)} rows and as many columns
## as there are inputs.  If @var{sys} is a single-input system, row vectors @var{u}
## of length @code{length(t)} are accepted as well.
## @item t
## Time vector.  Should be evenly spaced.  If @var{sys} is a continuous-time system
## and @var{t} is a real scalar, @var{sys} is discretized with sampling time
## @code{tsam = t/(rows(u)-1)}.  If @var{sys} is a discrete-time system and @var{t}
## is not specified, vector @var{t} is assumed to be @code{0 : tsam : tsam*(rows(u)-1)}.
## @item x0
## Vector of initial conditions for each state.  If not specified, a zero vector is assumed.
## @item method
## Discretization method for continuous-time models.  Default value is zoh
## (zero-order hold).  All methods from @code{c2d} are supported. 
## @end table
##
## @strong{Outputs}
## @table @var
## @item y
## Output response array.  Has as many rows as time samples (length of t)
## and as many columns as outputs.
## @item t
## Time row vector.  It is always evenly spaced.
## @item x
## State trajectories array.  Has @code{length (t)} rows and as many columns as states.
## @end table
##
## @seealso{impulse, initial, step}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.3

% function [y_r, t_r, x_r] = lsim (sys, u, t = [], x0 = [], method = "zoh")
function [y_r, t_r, x_r] = lsim (varargin)

  if (nargin < 2)
    print_usage ();
  endif

  sys_idx = cellfun (@isa, varargin, {"lti"});                          # look for LTI models
  sys_cell = cellfun (@ss, varargin(sys_idx), "uniformoutput", false);  # convert to state-space

  if (! size_equal (sys_cell{:}))
    error ("lsim: models must have equal sizes");
  endif

  vec_idx = find (cellfun (@is_real_matrix, varargin));
  n_vec = length (vec_idx);
  n_sys = length (sys_cell);



  if (! isa (sys, "ss"))
    sys = ss (sys);                             # sys must be proper
  endif

  if (is_real_vector (u))
    u = reshape (u, [], 1);                     # allow row vectors for single-input systems
  elseif (! is_real_matrix (u))
    error ("lsim: input signal u must be an array of real numbers");
  endif
  
  if (! is_real_vector (t) && ! isempty (t))
    error ("lsim: time vector t must be real or empty");
  endif
  
  discrete = ! isct (sys);                      # static gains are treated as continuous-time systems
  tsam = abs (get (sys, "tsam"));               # use 1 second as default if tsam is unspecified (-1)
  urows = rows (u);
  ucols = columns (u);

  if (discrete)                                 # discrete system
    if (isempty (t))                            # lsim (sys, u)
      dt = tsam;
      tinitial = 0;
      tfinal = tsam * (urows - 1);
    elseif (length (t) == 1)                    # lsim (sys, u, tfinal)
      dt = tsam;
      tinitial = 0;
      tfinal = t;
    else                                        # lsim (sys, u, t, ...)
      warning ("lsim: spacing of time vector has no effect on sampling time of discrete system");
      dt = tsam;
      tinitial = t(1);
      tfinal = t(end);
    endif
  else                                          # continuous system
    if (isempty (t))                            # lsim (sys, u, [], ...)
      error ("lsim: invalid time vector");
    elseif (length (t) == 1)                    # lsim (sys, u, tfinal, ...)
      dt = t / (urows - 1);
      tinitial = 0;
      tfinal = t;
    else                                        # lsim (sys, u, t, ...)
      dt = t(2) - t(1);                         # assume that t is regularly spaced
      tinitial = t(1);
      tfinal = t(end);
    endif
    sys = c2d (sys, dt, method);                # convert to discrete-time model
  endif

  [A, B, C, D] = ssdata (sys);

  n = rows (A);                                 # number of states
  m = columns (B);                              # number of inputs
  p = rows (C);                                 # number of outputs

  t = reshape (tinitial : dt : tfinal, [], 1);  # time vector
  trows = length (t);

  if (urows != trows)
    error ("lsim: input vector u must have %d rows", trows);
  endif
  
  if (ucols != m)
    error ("lsim: input vector u must have %d columns", m);
  endif

  ## preallocate memory
  y = zeros (trows, p);
  x_arr = zeros (trows, n);

  ## initial conditions
  if (isempty (x0))
    x0 = zeros (n, 1);
  elseif (n != length (x0) || ! is_real_vector (x0))
    error ("lsim: x0 must be a vector with %d elements", n);
  endif

  x = reshape (x0, [], 1);                      # make sure that x is a column vector

  ## simulation
  for k = 1 : trows
    y(k, :) = C * x  +  D * u(k, :).';
    x_arr(k, :) = x;
    x = A * x  +  B * u(k, :).';
  endfor

  if (nargout == 0)                             # plot information
    str = ["Linear Simulation Results of ", inputname(1)];
    outname = get (sys, "outname");
    outname = __labels__ (outname, "y_");
    if (discrete)                               # discrete system
      for k = 1 : p
        subplot (p, 1, k);
        stairs (t, y(:, k));
        grid ("on");
        if (k == 1)
          title (str);
        endif
        ylabel (sprintf ("Amplitude %s", outname{k}));
      endfor
      xlabel ("Time [s]");
    else                                        # continuous system
      for k = 1 : p
        subplot (p, 1, k);
        plot (t, y(:, k));
        grid ("on");
        if (k == 1)
          title (str);
        endif
        ylabel (sprintf ("Amplitude %s", outname{k}));
      endfor
      xlabel ("Time [s]");
    endif
  else                                          # return values
    y_r = y;
    t_r = t;
    x_r = x_arr;
  endif

endfunction


## TODO: add test cases