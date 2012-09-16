## Copyright (C) 2009, 2010, 2012   Lukas F. Reichlin
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
## Common code for the time response functions step, impulse and initial.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.3

% function [y, t, x_arr] = __time_response_2__ (sys, resptype, plotflag, tfinal, dt, x0, sysname)
function [y, t, x] = __time_response_2__ (resptype, args)

  sys_idx = cellfun (@isa, args, {"lti"});      # look for LTI models
  sys_cell = cellfun (@ss, args(sys_idx));      # system must be proper

  if (! size_equal (sys_cell{:}))
    error ("%s: models must have equal sizes", resptype);
  endif

  tmp = cellfun (@is_real_matrix, args);
  vec_idx = find (tmp);
  n_vec = length (vec_idx)
  n_sys = length (sys_cell)
  
  if (n_vec >= 1)
    arg = args{vec_idx(1)};
    
  endif

  ## extract tfinal/t, dt, x0

  [tfinal, dt] = cellfun (@__sim_horizon__, sys_cell, {tfinal}, dt);
  
  tfinal = max (tfinal);
 
  % __sim_horizon__ (sys, tfinal, dt);
 
 
  hier sim_horizon
  
  ct_idx = cellfun (@isct, sys_cell)
  sys_dt_cell = sys_cell;
  tmp = (@c2d, sys_cell(ct_idx), dt, {"zoh"}, "uniformoutput", false)
  sys_dt_cell(ct_idx) = tmp;

%{
  if (! isa (sys, "ss"))
    sys = ss (sys);                                    # sys must be proper
  endif

  if (is_real_vector (tfinal) && length (tfinal) > 1)  # time vector t passed
    dt = tfinal(2) - tfinal(1);                        # assume that t is regularly spaced
    tfinal = tfinal(end);
  endif

  [A, B, C, D, tsam] = ssdata (sys);

  discrete = ! isct (sys);                             # static gains are treated as analog systems
  tsam = abs (tsam);                                   # use 1 second if tsam is unspecified (-1)

  if (discrete)
    if (! isempty (dt))
      warning ("time_response: argument dt has no effect on sampling time of discrete system");
    endif

    dt = tsam;
  endif

  [tfinal, dt] = __sim_horizon__ (A, discrete, tfinal, dt);

  if (! discrete)
    sys = c2d (sys, dt, "zoh");
  endif

  [F, G] = ssdata (sys);                               # matrices C and D don't change

  n = rows (F);                                        # number of states
  m = columns (G);                                     # number of inputs
  p = rows (C);                                        # number of outputs

  ## time vector
  t = reshape (0 : dt : tfinal, [], 1);
  l_t = length (t);
%}

%function [y, x_arr] = __initial_response__ (sys, sys_dt, t, x0)
%function [y, x_arr] = __step_response__ (sys_dt, t)
%function [y, x_arr] = __impulse_response__ (sys, sys_dt, t)

  switch (resptype)
    case "initial"
      [y, x] = cellfun (@__initial_response__, sys_cell, sys_dt_cell, t, {x0} or x0);
    case "step"
      [y, x] = cellfun (@__step_response__, sys_dt_cell, t);
    case "impulse"
      [y, x] = cellfun (@__impulse_response__, sys_cell, sys_dt_cell, t);
    otherwise
      error ("time_response: invalid response type");
  endswitch




  if (plotflag)                                        # display plot
    switch (resptype)
      case "initial"
        str = "Response to Initial Conditions";
        cols = 1;
      case "step"
        str = "Step Response";
        cols = m;
      case "impulse"
        str = "Impulse Response";
        cols = m;
      otherwise
        error ("time_response: invalid response type");
    endswitch
  
    for i = 1 : n_sys
      t = t{i};
      y = y{i};
      discrete = ! ct_idx(i);
      if (discrete)                                      # discrete system
        for k = 1 : p
          for j = 1 : cols
            subplot (p, cols, (k-1)*cols+j);
            stairs (t, y(:, k, j));
            hold on;
            grid ("on");
            if (i == n_sys && k == 1 && j == 1)
              title (str);
            endif
            if (i == n_sys && j == 1)
              ylabel (sprintf ("Amplitude %s", outname{k}));
            endif
          endfor
        endfor
      else                                               # continuous system
        for k = 1 : p
          for j = 1 : cols
            subplot (p, cols, (k-1)*cols+j);
            plot (t, y(:, k, j));
            hold on;
            grid ("on");
            if (i == n_sys && k == 1 && j == 1)
              title (str);
            endif
            if (i == n_sys && j == 1)
              ylabel (sprintf ("Amplitude %s", outname{k}));
            endif
          endfor
        endfor
      endif
    endfor
    xlabel ("Time [s]");
    hold off;
  endif

%{  
  if (plotflag)                                        # display plot

    ## TODO: Set correct titles, especially for multi-input systems

    stable = isstable (sys);
    outname = get (sys, "outname");
    outname = __labels__ (outname, "y_");

    if (strcmp (resptype, "initial"))
      cols = 1;
    else
      cols = m;
    endif

    if (discrete)                                      # discrete system
      for k = 1 : p
        for j = 1 : cols

          subplot (p, cols, (k-1)*cols+j);

          if (stable)
            stairs (t, [y(:, k, j), yfinal(k, j) * ones(l_t, 1)]);
          else
            stairs (t, y(:, k, j));
          endif

          grid ("on");

          if (k == 1 && j == 1)
            title (str);
          endif

          if (j == 1)
            ylabel (sprintf ("Amplitude %s", outname{k}));
          endif

        endfor
      endfor

      xlabel ("Time [s]");

    else                                               # continuous system
      for k = 1 : p
        for j = 1 : cols

          subplot (p, cols, (k-1)*cols+j);

          if (stable)
            plot (t, [y(:, k, j), yfinal(k, j) * ones(l_t, 1)]);
          else
            plot (t, y(:, k, j));
          endif

          grid ("on");

          if (k == 1 && j == 1)
            title (str);
          endif

          if (j == 1)
            ylabel (sprintf ("Amplitude %s", outname{k}));
          endif

        endfor
      endfor

      xlabel ("Time [s]");

    endif 
  endif

endfunction
%}

endfunction




function [y, x_arr] = __initial_response__ (sys, sys_dt, t, x0)

  [A, B, C, D] = ssdata (sys);
  [F, G] = ssdata (sys_dt);

  n = rows (F);                                        # number of states
  m = columns (G);                                     # number of inputs
  p = rows (C);                                        # number of outputs
  l_t = length (t);

  ## preallocate memory
  y = zeros (l_t, p);
  x_arr = zeros (l_t, n);

  ## initial conditions
  x = reshape (x0, [], 1);                         # make sure that x is a column vector

  if (n != length (x0) || ! is_real_vector (x0))
    error ("initial: x0 must be a real vector with %d elements", n);
  endif

  ## simulation
  for k = 1 : l_t
    y(k, :) = C * x;
    x_arr(k, :) = x;
    x = F * x;
  endfor

endfunction
  

function [y, x_arr] = __step_response__ (sys_dt, t)

  [F, G, C, D] = ssdata (sys_dt);

  n = rows (F);                                        # number of states
  m = columns (G);                                     # number of inputs
  p = rows (C);                                        # number of outputs
  l_t = length (t);

  ## preallocate memory
  y = zeros (l_t, p, m);
  x_arr = zeros (l_t, n, m);

  for j = 1 : m                                    # for every input channel
    ## initial conditions
    x = zeros (n, 1);
    u = zeros (m, 1);
    u(j) = 1;

    ## simulation
    for k = 1 : l_t
      y(k, :, j) = C * x + D * u;
      x_arr(k, :, j) = x;
      x = F * x + G * u;
    endfor
  endfor

endfunction


function [y, x_arr] = __impulse_response__ (sys, sys_dt, t)

  [A, B, C, D, dt] = ssdata (sys);
  [F, G] = ssdata (sys_dt);
  discrete = ! isct (sys);

  n = rows (F);                                        # number of states
  m = columns (G);                                     # number of inputs
  p = rows (C);                                        # number of outputs
  l_t = length (t);

  ## preallocate memory
  y = zeros (l_t, p, m);
  x_arr = zeros (l_t, n, m);

  for j = 1 : m                                    # for every input channel
    ## initial conditions
    u = zeros (m, 1);
    u(j) = 1;

    if (discrete)
      x = zeros (n, 1);                            # zero by definition 
      y(1, :, j) = D * u / dt;
      x_arr(1, :, j) = x;
      x = G * u / dt;
    else
      x = B * u;                                   # B, not G!
      y(1, :, j) = C * x;
      x_arr(1, :, j) = x;
      x = F * x;
    endif

    ## simulation
    for k = 2 : l_t
      y (k, :, j) = C * x;
      x_arr(k, :, j) = x;
      x = F * x;
    endfor
  endfor

  if (discrete)
    y *= dt;
    x_arr *= dt;
  endif

endfunction


function 



% function [tfinal, dt] = __sim_horizon__ (A, discrete, tfinal, Ts)
function [tfinal, dt] = __sim_horizon__ (sys, tfinal, Ts)

  ## code based on __stepimp__.m of Kai P. Mueller and A. Scottedward Hodel

  TOL = 1.0e-10;                                       # values below TOL are assumed to be zero
  N_MIN = 50;                                          # min number of points
  N_MAX = 2000;                                        # max number of points
  N_DEF = 1000;                                        # default number of points
  T_DEF = 10;                                          # default simulation time

  ev = pole (sys);
  n = length (ev);

  if (discrete)
    ## perform bilinear transformation on poles in z
    for k = 1 : n
      pol = ev(k);
      if (abs (pol + 1) < TOL)
        ev(k) = 0;
      else
        ev(k) = 2 / Ts * (pol - 1) / (pol + 1);
      endif
    endfor
  endif

  ## remove poles near zero from eigenvalue array ev
  nk = n;
  for k = 1 : n
    if (abs (real (ev(k))) < TOL)
      ev(k) = 0;
      nk -= 1;
    endif
  endfor

  if (nk == 0)
    if (isempty (tfinal))
      tfinal = T_DEF;
    endif

    if (! discrete)
      dt = tfinal / N_DEF;
    endif
  else
    ev = ev(find (ev));
    ev_max = max (abs (ev));

    if (! discrete)
      dt = 0.2 * pi / ev_max;
    endif

    if (isempty (tfinal))
      ev_min = min (abs (real (ev)));
      tfinal = 5.0 / ev_min;

      ## round up
      yy = 10^(ceil (log10 (tfinal)) - 1);
      tfinal = yy * ceil (tfinal / yy);
    endif

    if (! discrete)
      N = tfinal / dt;

      if (N < N_MIN)
        dt = tfinal / N_MIN;
      endif

      if (N > N_MAX)
        dt = tfinal / N_MAX;
      endif
    endif
  endif

  if (! isempty (Ts))                                  # catch case cont. system with dt specified
    dt = Ts;
  endif

endfunction
