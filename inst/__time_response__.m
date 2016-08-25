## Copyright (C) 2009-2016   Lukas F. Reichlin
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
## Version: 0.5

function [y, t, x] = __time_response__ (response, args, nout)

  idx = cellfun (@islogical, args);
  tmp = cellfun (@double, args(idx), "uniformoutput", false);
  args(idx) = tmp;

  sys_idx = cellfun (@isa, args, {"lti"});                          # LTI models
  mat_idx = cellfun (@is_real_matrix, args);                        # matrices
  sty_idx = cellfun (@ischar, args);                                # strings (style arguments)

  inv_idx = ! (sys_idx | mat_idx | sty_idx);                        # invalid arguments

  if (any (inv_idx))
    warning ("%s: arguments number %s are invalid and are being ignored", ...
             response, mat2str (find (inv_idx)(:).'));
  endif

  if (nnz (sys_idx) == 0)
    error ("%s: require at least one LTI model", response);
  endif

  if (nout > 0 && (nnz (sys_idx) > 1 || any (sty_idx)))
    evalin ("caller", "print_usage ()");
  endif

  if (! size_equal (args{sys_idx}))
    error ("%s: all LTI models must have equal size", response);
  endif

  if (any (find (sty_idx) < find (sys_idx)(1)))
    warning ("%s: strings in front of first LTI model are being ignored", response);
  endif

  tfinal = [];  dt = [];  x0 = [];                                  # default arguments

  switch (response)
    case "initial"
      switch (nnz (mat_idx))
        case 0
          error ("initial: require initial state vector 'x0'");
        case 1
          x0 = args{mat_idx};
        case 2
          [x0, tfinal] = args{mat_idx};
        case 3
          [x0, tfinal, dt] = args{mat_idx};
        otherwise
          evalin ("caller", "print_usage ()");
      endswitch
      if (! is_real_vector (x0))
        error ("initial: initial state vector 'x0' must be a real-valued vector");
      endif
    
    case {"step", "impulse", "ramp"}
      switch (nnz (mat_idx))
        case 0
          ## nothing to here, just prevent case 'otherwise'
        case 1
          tfinal = args{mat_idx};
        case 2
          [tfinal, dt] = args{mat_idx};
        otherwise
          evalin ("caller", "print_usage ()");
      endswitch
    
    otherwise
      error ("time_response: invalid response type '%s'", response);
  endswitch
  
  switch (response)
    case "step"
      response1 = "zoh";
    case "impulse"
      response1 = "impulse";
    otherwise
      response1 = "zoh";
  endswitch
  
  if (issample (tfinal) || isempty (tfinal))
    ## nothing to do here
  elseif (is_real_vector (tfinal))
    dt = abs (tfinal(end) - tfinal(1)) / (length (tfinal) - 1);
    tfinal = abs (tfinal(end));
  else
    evalin ("caller", "print_usage ()");
  endif

  if (isempty (dt))
    ## nothing to do here
  elseif (issample (dt))
    ## nothing to do here
  else
    evalin ("caller", "print_usage ()");
  endif

  [tfinal, dt] = cellfun (@__sim_horizon__, args(sys_idx), {tfinal}, {dt}, "uniformoutput", false);
  tfinal = max ([tfinal{:}]);

  ct_idx = cellfun (@isct, args(sys_idx));
  sys_dt = args(sys_idx);
  tmp = cellfun (@c2d, args(sys_idx)(ct_idx), dt(ct_idx), {response1}, "uniformoutput", false);
  sys_dt(ct_idx) = tmp;

  ## time vector
  t = cellfun (@(dt) vec (linspace (0, tfinal, tfinal/dt+1)), dt, "uniformoutput", false);
  
  ## alternative code
  ## t = cellfun (@(dt) vec (0 : dt : tfinal), dt, "uniformoutput", false);

  ## function [y, x_arr] = __initial_response__ (sys_dt, t, x0)
  ## function [y, x_arr] = __step_response__ (sys_dt, t)
  ## function [y, x_arr] = __impulse_response__ (sys, sys_dt, t)
  ## function [y, x_arr] = __ramp_response__ (sys_dt, t)

  switch (response)
    case "initial"
      [y, x] = cellfun (@__initial_response__, sys_dt, t, {x0}, "uniformoutput", false);
    case "step"
      [y, x] = cellfun (@__step_response__, sys_dt, t, "uniformoutput", false);
    case "impulse"
      [y, x] = cellfun (@__impulse_response__, args(sys_idx), sys_dt, t, "uniformoutput", false);
    case "ramp"
      [y, x] = cellfun (@__ramp_response__, sys_dt, t, "uniformoutput", false);
    otherwise
      error ("time_response: invalid response type");
  endswitch


  if (nout == 0)                                        # display plot
    ## extract plotting styles
    tmp = cumsum (sys_idx);
    tmp(sys_idx | ! sty_idx) = 0;
    n_sys = nnz (sys_idx);
    sty = arrayfun (@(x) args(tmp == x), 1:n_sys, "uniformoutput", false);

    ## default plotting styles if empty
    colororder = get (gca, "colororder");
    rc = rows (colororder);
    def = arrayfun (@(k) {"color", colororder(1+rem (k-1, rc), :)}, 1:n_sys, "uniformoutput", false);
    idx = cellfun (@isempty, sty);
    sty(idx) = def(idx);

    ## get system names for legend
    leg = cell (1, n_sys);
    idx = find (sys_idx);
    for k = 1 : n_sys
      leg{k} = evalin ("caller", sprintf ("inputname(%d)", idx(k)), "''");
    endfor

    outname = get (args(sys_idx){end}, "outname");
    outname = __labels__ (outname, "y");

    [p, m] = size (args(sys_idx){1});
    
    switch (response)
      case "initial"
        str = "Response to Initial Conditions";
        cols = 1;
        ## yfinal = zeros (p, 1);
      case "step"
        str = "Step Response";
        cols = m;
        ## yfinal = dcgain (sys_cell{1});
      case "impulse"
        str = "Impulse Response";
        cols = m;
        ## yfinal = zeros (p, m);
      case "ramp"
        str = "Ramp Response";
        cols = m;
      otherwise
        error ("time_response: invalid response type");
    endswitch


    for k = 1 : n_sys                                   # for every system
      if (ct_idx(k))                                    # continuous-time system                                           
        for i = 1 : p                                   # for every output
          for j = 1 : cols                              # for every input (except for initial where cols=1)
            if (p != 1 || cols != 1)
              subplot (p, cols, (i-1)*cols+j);
            endif
            plot (t{k}, y{k}(:, i, j), sty{k}{:});
            hold on;
            grid on;
            if (k == n_sys)
              axis tight
              ylim (__axis_margin__ (ylim))
              if (j == 1)
                ylabel (outname{i});
                if (i == 1)
                  title (str);
                endif
              endif
            endif
          endfor
        endfor
      else                                              # discrete-time system
        for i = 1 : p                                   # for every output
          for j = 1 : cols                              # for every input (except for initial where cols=1)
            if (p != 1 || cols != 1)
              subplot (p, cols, (i-1)*cols+j);
            endif
            stairs (t{k}, y{k}(:, i, j), sty{k}{:});
            hold on;
            grid on;
            if (k == n_sys)
              axis tight;
              ylim (__axis_margin__ (ylim))
              if (j == 1)
                ylabel (outname{i});
                if (i == 1)
                  title (str);
                endif
              endif
            endif
          endfor
        endfor
      endif
    endfor
    xlabel ("Time [s]");
    if (p == 1 && m == 1)
      legend (leg)
    endif
    hold off;
  endif

endfunction


function [y, x_arr] = __initial_response__ (sys_dt, t, x0)

  [F, G, C, D] = ssdata (sys_dt);                       # system must be proper

  n = rows (F);                                         # number of states
  m = columns (G);                                      # number of inputs
  p = rows (C);                                         # number of outputs
  l_t = length (t);

  ## preallocate memory
  y = zeros (l_t, p);
  x_arr = zeros (l_t, n);

  ## initial conditions
  x = reshape (x0, [], 1);                              # make sure that x is a column vector

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

  [F, G, C, D] = ssdata (sys_dt);       # system must be proper

  n = rows (F);                                         # number of states
  m = columns (G);                                      # number of inputs
  p = rows (C);                                         # number of outputs
  l_t = length (t);

  ## preallocate memory
  y = zeros (l_t, p, m);
  x_arr = zeros (l_t, n, m);

  for j = 1 : m                                         # for every input channel
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

 # [~, B] = ssdata (sys);
  [F, G, C, D, dt] = ssdata (sys_dt);                   # system must be proper
  dt = abs (dt);                                        # use 1 second if tsam is unspecified (-1)
  discrete = ! isct (sys_dt);

  n = rows (F);                                         # number of states
  m = columns (G);                                      # number of inputs
  p = rows (C);                                         # number of outputs
  l_t = length (t);

  ## preallocate memory
  y = zeros (l_t, p, m);
  x_arr = zeros (l_t, n, m);

  for j = 1 : m                                         # for every input channel
    ## initial conditions
    u = zeros (m, 1);
    u(j) = 1;

    if (discrete)
      x = zeros (n, 1);                                 # zero by definition 
      y(1, :, j) = D * u / dt;
      x_arr(1, :, j) = x;
      x = G * u / dt;
    else
      x = G * u;                                        #NO NO B, not G!
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


function [y, x_arr] = __ramp_response__ (sys_dt, t)

  [F, G, C, D] = ssdata (sys_dt);       # system must be proper

  n = rows (F);                                         # number of states
  m = columns (G);                                      # number of inputs
  p = rows (C);                                         # number of outputs
  l_t = length (t);

  ## preallocate memory
  y = zeros (l_t, p, m);
  x_arr = zeros (l_t, n, m);

  for j = 1 : m                                         # for every input channel
    ## initial conditions
    x = zeros (n, 1);
    u = zeros (m, l_t);
    u(j, :) = t;

    ## simulation
    for k = 1 : l_t
      y(k, :, j) = C * x + D * u(:, k);
      x_arr(k, :, j) = x;
      x = F * x + G * u(:, k);
    endfor
  endfor

endfunction


function [tfinal, dt] = __sim_horizon__ (sys, tfinal, Ts)

  ## code based on __stepimp__.m of Kai P. Mueller and A. Scottedward Hodel

  TOL = 1.0e-10;                                        # values below TOL are assumed to be zero
  N_MIN = 50;                                           # min number of points
  N_MAX = 2000;                                         # max number of points
  N_DEF = 1000;                                         # default number of points
  T_DEF = 10;                                           # default simulation time

  ev = pole (sys);
  n = length (ev);                                      # number of states/poles
  continuous = isct (sys);
  discrete = ! continuous;

  if (discrete)
    dt = Ts = abs (get (sys, "tsam"));
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

    if (continuous)
      dt = tfinal / N_DEF;
    endif
  else
    ev = ev(find (ev));
    ev_max = max (abs (ev));

    if (continuous)
      dt = 0.2 * pi / ev_max;
    endif

    if (isempty (tfinal))
      ev_min = min (abs (real (ev)));
      tfinal = 5.0 / ev_min;

      ## round up
      yy = 10^(ceil (log10 (tfinal)) - 1);
      tfinal = yy * ceil (tfinal / yy);
    endif

    if (continuous)
      N = tfinal / dt;

      if (N < N_MIN)
        dt = tfinal / N_MIN;
      endif

      if (N > N_MAX)
        dt = tfinal / N_MAX;
      endif
    endif
  endif

  if (continuous && ! isempty (Ts))                     # catch case cont. system with dt specified
    dt = Ts;
  endif

endfunction
