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

function [y, t, x] = __time_response__ (response, args, nout=0, names={})
  ## varargin: nout = 0, names

  idx = cellfun (@islogical, args);
  tmp = cellfun (@double, args(idx), "uniformoutput", false);
  args(idx) = tmp;

  [sys_idx, sty_idx, mat_idx, names, sty] = ...
      __control_args__ (args, {"@lti", @ischar, @is_real_matrix}, names);

  inv_idx = ! (sys_idx | mat_idx | sty_idx);      # invalid arguments

  if (any (inv_idx))
    warning ("%s: arguments number %s are invalid and are being ignored\n", ...
             response, mat2str (find (inv_idx)(:).'));
  endif

  n_sys = nnz (sys_idx);
  if (n_sys == 0)
    error ("%s: require at least one LTI model\n", response);
  endif

  if (nout > 0)
    if nnz (sys_idx) > 1
      error ("%s: with output arguments, only one system is allowed\n", response);
    endif
    if any (sty_idx)
      warning ("%s: with output arguments, all style parameters are ignored\n", response);
    endif
  endif

  if (any (find (sty_idx) < find (sys_idx)(1)))
    warning ("%s: strings in front of first LTI model are being ignored\n", response);
  endif

    tfinal = [];  dt = [];  x0 = [];                                  # default arguments

  switch (response)
    case "initial"
      switch (nnz (mat_idx))
        case 0
          error ("initial: require initial state vector 'x0'\n");
        case 1
          x0 = args{mat_idx};
        case 2
          [x0, tfinal] = args{mat_idx};
        case 3
          [x0, tfinal, dt] = args{mat_idx};
        otherwise
          print_usage (response);
      endswitch
      if (! is_real_vector (x0))
        error ("initial: initial state vector 'x0' must be a real-valued vector\n");
      endif
      is_ss = cellfun (@isa, args(sys_idx), {"ss"});
      if (! all (is_ss))
        no_ss = find (is_ss==0);
        no_ss_names = sprintf ("%s, ", names{no_ss});
        no_ss_names = no_ss_names(1:end-2);
        error ("initial: system %s not in state space, x0 is ambiguous\n", no_ss_names);
      endif

    case {"step", "stepinfo", "impulse", "ramp"}
      switch (nnz (mat_idx))
        case 0
          ## nothing to here, just prevent case 'otherwise'
        case 1
          tfinal = args{mat_idx};
        case 2
          [tfinal, dt] = args{mat_idx};
        otherwise
          print_usage (response);
      endswitch

    otherwise
      error ("time_response: invalid response type '%s'\n", response);
  endswitch

  switch (response)
    case {"step", "stepresponse"}
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
    print_usage (response);
  endif

  if (isempty (dt))
    ## nothing to do here
  elseif (issample (dt))
    ## nothing to do here
  else
    print_usage (response);
  endif

  [tfinal, dt] = cellfun (@__sim_horizon__, args(sys_idx), {tfinal}, {dt}, "uniformoutput", false);
  tfinal = max ([tfinal{:}]);

  ## discretizaiton of continuous time systems
  ## do this in state space for more accurate results
  sys_dt = args(sys_idx);
  ct_idx = cellfun (@isct, sys_dt);
  sys_ct = sys_dt(ct_idx);
  sys_idx_ct = find (ct_idx);
  ## FIXME: ss can not be applied via cellfun ()? Use a for-loop instead
  ##        "lti: subsasgn: invalid subscripted assignment type '()'"
  sys_ctss = cell (size (sys_ct));
  for i = 1:length (sys_ct)
    sys_ctss{i} = ss (sys_ct{i});
    # impulse response for systems with direct feedthrough needs special care:
    # get the response without the resulting delta-inpulses in the outputs
    # and print a warning
    if strcmp (response1,"impulse")
      [nz_y, nz_u] = find (sys_ctss{i}.d);
      if ! isempty (nz_y)
        # there is at least one direct feedthrough, prepare a warning message
        msg = ["System \"%s\" has direct feedthrough!\n", ...
               "The impulse %f*delta(t) is omitted in the impulse response ", ...
               "from input \"%s\" to output \"%s\"\n"];
      endif
      for jy = 1:length (nz_y)
        for ju = 1:length (nz_u)
          # get the names of all input/output channels with feedthrough
          in = sys_ctss{i}.inputname{nz_u(ju)};
          out = sys_ctss{i}.outputname{nz_y(jy)};
          if isempty (in)
            in = ['u', num2str(nz_u(ju))];  # default name
          endif
          if isempty (out)
            out = ['y', num2str(nz_y(jy))];  # default name
          endif
          # print the warning
          warning (msg, names{sys_idx_ct(i)}, sys_ctss{i}.d(nz_y(jy),nz_u(ju)), in, out);
        endfor
      endfor
      ## compute impulse response for remaining system without feedthrough
      sys_ctss{i}.d = 0*sys_ctss{i}.d;
    endif
  endfor
  sys_ct2dt = cellfun (@c2d, sys_ctss, dt(ct_idx), {response1}, "uniformoutput", false);
  sys_dt(ct_idx) = sys_ct2dt;

  ## delay is already whole samples on sys_dt (c2d's default DelayModeling
  ## rounds to samples at conversion time; already-discrete inputs are
  ## documented in samples too), so every shift below is an exact
  ## integer-sample array shift -- no interpolation needed.
  delay_mat = cellfun (@totaldelay, sys_dt, "uniformoutput", false);

  ## time vector: we have to consider the following cases:
  ##              1. ct system: last sample is tfinal (ensured by __sim_horizon__)
  ##              2. dt system
  ##                  a) nout > 0 (no plotting): last sample is less or equal tfinal
  ##                  b) nout > 0 (plotting): last sample is the first greater
  ##                     than tfinal (we need xlim([0,tfinal]) for the plot)
  delay_ext = 0;
  for i = 1:length (dt)
    delay_ext = max (delay_ext, max (delay_mat{i}(:)) * dt{i});
  endfor
  tfinal = tfinal + delay_ext;

  if nout > 0
    dt_extra = cell2mat (dt) .* ct_idx;
  else
    dt_extra = cell2mat (dt);
  end
  factor_horizon = 1;
  if (strcmp (response, "stepinfo"))
    factor_horizon = 1.5;   # longer horizon for stepinfo
  endif
  t = cell (size(dt));
  for i = 1:length(t)
    t{i} = vec (0:dt{i}:factor_horizon*tfinal);
    if (ct_idx(i) == 0) && (nout == 0) && (length (t{i}) * dt{i} < tfinal)
      ## Discrete time system, no plotting, and last sampling is before tfinal
      t{i}(end+1) = t{i}(end) + dt{i};
    end
  end

  ## alternative code
  ## t = cellfun (@(dt) vec (0 : dt : tfinal), dt, "uniformoutput", false);

  ## function [y, x_arr] = __initial_response__ (sys_dt, t, x0)
  ## function [y, x_arr] = __step_response__ (sys_dt, t)
  ## function [y, x_arr] = __impulse_response__ (sys, sys_dt, t)
  ## function [y, x_arr] = __ramp_response__ (sys_dt, t)

  switch (response)
    case "initial"
      [y, x] = cellfun (@__initial_response__, sys_dt, t, {x0}, "uniformoutput", false);
    case {"step", "stepinfo"}
      [y, x] = cellfun (@__step_response__, sys_dt, t, "uniformoutput", false);
    case "impulse"
      [y, x] = cellfun (@__impulse_response__, args(sys_idx), sys_dt, t, "uniformoutput", false);
    case "ramp"
      [y, x] = cellfun (@__ramp_response__, sys_dt, t, "uniformoutput", false);
    otherwise
      error ("time_response: invalid response type\n");
  endswitch

  ## apply each channel's total delay as a post-processing shift on y;
  ## x (state trajectory) is intentionally left unshifted, see plan/spec.
  for i = 1:n_sys
    if (strcmp (response, "initial"))
      outdelay = get (sys_dt{i}, "outputdelay");
      for row = 1:columns (y{i})
        y{i}(:, row) = __apply_timeresp_delay__ (y{i}(:, row), outdelay(row));
      endfor
    else
      p_i = size (y{i}, 2);
      m_i = size (y{i}, 3);
      for row = 1:p_i
        for col = 1:m_i
          y{i}(:, row, col) = __apply_timeresp_delay__ (y{i}(:, row, col), delay_mat{i}(row, col));
        endfor
      endfor
    endif
  endfor


  if (nout == 0)                                        # display plot

    ## get max sizes
    [p, m] = cellfun (@size, args(sys_idx), 'uniformoutput', false);
    p = cell2mat (p);
    m = cell2mat (m);
    rows = max (p);
    cols = max (m);

    ## get in/outnames, take ui and yi for more than one system
    if n_sys == 1
      outname = get (args(sys_idx){1}, "outname");
      inname = get (args(sys_idx){1}, "inname");
    else
      outname = cell(1,rows);
      inname = cell(1,cols);
    endif
    outname = __labels__ (outname, "y");
    inname = __labels__ (inname, "u");

    switch (response)
      case "initial"
        str = "Response to Initial Conditions";
        cols = 1;
        ## yfinal = zeros (p, 1);
      case {"step", "stepinfo"}
        str = "Step Response";
        ## yfinal = dcgain (sys_cell{1});
      case "impulse"
        str = "Impulse Response";
        ## yfinal = zeros (p, m);
      case "ramp"
        str = "Ramp Response";
      otherwise
        error ("time_response: invalid response type '%s'\n", response);
    endswitch

    ## get last system present in the subplots
    last_system = zeros (rows, cols);
    for i = 1 : rows
      for j = 1 : cols
        for k = 1 : n_sys
          if (p(k) >= i) && (m(k) >= j)
            last_system(i,j) = k;
          endif
        endfor
      endfor
    endfor

    ## loop over all subplots
    for i = 1 : rows              # for every output
      for j = 1 : cols            # for every input (except for initial where cols=1)
        if last_system(i,j) > 0   # only if there is a system with this in-/output combination

          if (rows != 1) || (cols != 1)
            ## if we only have one plot, don't use subplot
            ## allowing a step response in an existing subplot
            subplot (rows, cols, (i-1)*cols+j);
          endif

          box on;

          for k = 1 : last_system(i,j)                # for every system for this in-/output

            p_sys = size (args(sys_idx){k}, 1);
            m_sys = size (args(sys_idx){k}, 2);
            if (p_sys >= i) && (m_sys >= j)           # only if system has this in-/output

                if (ct_idx(k))    # continuous-time system
                  plot (t{k}, y{k}(:, i, j), sty{k}{:});
                else              # discrete-time system
                  [tstep,ystep] = stairs (t{k}, y{k}(:, i, j));
                  plot (tstep, ystep, sty{k}{:});
                endif
                hold on;
                if k == last_system(i,j)
                  ## do the following only for last system in this subplot
                  grid on;
                  axis tight
                  xlim ([0, tfinal]);
                  ylim (__axis_margin__ (ylim))
                  xlabel ("Time [s]");
                  ylabel (outname{i});
                  if (! strcmp (response, "initial")) && (rows > 1 || cols > 1)
                    title (inname{j}, "fontweight", "normal");
                  endif
                endif

            endif

          endfor

          hold off;

        endif

        hold off;

      endfor

    endfor

    if (rows == 1 && cols == 1)
      title (str);  # normal title
    else
      # create title manually (create axes object on whole figure and put text)
      f_pos = get (gcf(), 'position');
      fs = get (gca(), 'fontsize');
      last_plot = gca ();
      dummy = axes( 'visible', 'off', 'position', [0 0 1 1]);
      text (dummy, 0.5, 1-1.2*fs/f_pos(4), str, ...
            'fontsize', 1.2*fs, 'fontweight', 'bold', 'horizontalalignment', 'center');
      set (gcf (), 'currentaxes', last_plot); # focus back to last subplot
    endif

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
    error ("initial: x0 must be a real vector with %d elements\n", n);
  endif

  if (hasinternaldelay (sys_dt))
    [B2, C2, D12, D21, D22, tau] = __internaldelay_ports__ (sys_dt);
    nports = numel (tau);
    z_hist = zeros (l_t, nports);
    ## simulation (input is zero; only the internal-delay port drives w/z)
    for k = 1 : l_t
      w = __delay_lookup__ (z_hist, k, tau, nports);
      y(k, :) = C * x + D12 * w;
      z_hist(k, :) = (C2 * x + D22 * w).';
      x_arr(k, :) = x;
      x = F * x + B2 * w;
    endfor
    return;
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

  has_id = hasinternaldelay (sys_dt);
  if (has_id)
    [B2, C2, D12, D21, D22, tau] = __internaldelay_ports__ (sys_dt);
    nports = numel (tau);
  endif

  for j = 1 : m                                         # for every input channel
    ## initial conditions
    x = zeros (n, 1);
    u = zeros (m, 1);
    u(j) = 1;

    if (has_id)
      z_hist = zeros (l_t, nports);
      for k = 1 : l_t
        w = __delay_lookup__ (z_hist, k, tau, nports);
        y(k, :, j) = C * x + D * u + D12 * w;
        z_hist(k, :) = (C2 * x + D21 * u + D22 * w).';
        x_arr(k, :, j) = x;
        x = F * x + G * u + B2 * w;
      endfor
    else
      ## simulation
      for k = 1 : l_t
        y(k, :, j) = C * x + D * u;
        x_arr(k, :, j) = x;
        x = F * x + G * u;
      endfor
    endif
  endfor

endfunction


function [y, x_arr] = __impulse_response__ (sys, sys_dt, t)

  [F, G, C, D, dt] = ssdata (sys_dt);                   # system must be proper
  dt = abs (dt);                                        # use 1 second if tsam is unspecified (-1)

  n = rows (F);                                         # number of states
  m = columns (G);                                      # number of inputs
  p = rows (C);                                         # number of outputs
  l_t = length (t);

  ## preallocate memory
  y = zeros (l_t, p, m);
  x_arr = zeros (l_t, n, m);

  has_id = hasinternaldelay (sys_dt);
  if (has_id)
    [B2, C2, D12, D21, D22, tau] = __internaldelay_ports__ (sys_dt);
    nports = numel (tau);
  endif

  for j = 1 : m                                         # for every input channel

    u = zeros (m, 1);
    u(j) = 1;

    if (has_id)
      ## Equivalent to the delay-free branch below: a discrete impulse is the
      ## input signal uu(1) = u/dt, uu(k>1) = 0, run through the standard
      ## delay-buffer recursion (the k=1 special-case there is just this loop
      ## evaluated at the first sample).
      x = zeros (n, 1);
      z_hist = zeros (l_t, nports);
      for k = 1 : l_t
        if (k == 1), uk = u / dt; else, uk = zeros (m, 1); endif
        w = __delay_lookup__ (z_hist, k, tau, nports);
        y(k, :, j) = C * x + D * uk + D12 * w;
        z_hist(k, :) = (C2 * x + D21 * uk + D22 * w).';
        x_arr(k, :, j) = x;
        x = F * x + G * uk + B2 * w;
      endfor
    else
      ## initial conditions
      x = zeros (n, 1);                                 # zero by definition
      y(1, :, j) = D * u / dt;                          # impulse is 1/dt
      x_arr(1, :, j) = x;
      x = G * u / dt;

      ## simulation
      for k = 2 : l_t
        y (k, :, j) = C * x;
        x_arr(k, :, j) = x;
        x = F * x;
      endfor
    endif

  endfor

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

  has_id = hasinternaldelay (sys_dt);
  if (has_id)
    [B2, C2, D12, D21, D22, tau] = __internaldelay_ports__ (sys_dt);
    nports = numel (tau);
  endif

  for j = 1 : m                                         # for every input channel
    ## initial conditions
    x = zeros (n, 1);
    u = zeros (m, l_t);
    u(j, :) = t;

    if (has_id)
      z_hist = zeros (l_t, nports);
      for k = 1 : l_t
        w = __delay_lookup__ (z_hist, k, tau, nports);
        y(k, :, j) = C * x + D * u(:, k) + D12 * w;
        z_hist(k, :) = (C2 * x + D21 * u(:, k) + D22 * w).';
        x_arr(k, :, j) = x;
        x = F * x + G * u(:, k) + B2 * w;
      endfor
    else
      ## simulation
      for k = 1 : l_t
        y(k, :, j) = C * x + D * u(:, k);
        x_arr(k, :, j) = x;
        x = F * x + G * u(:, k);
      endfor
    endif
  endfor

endfunction


## Extract the internal-delay port matrices (B2, C2, D12, D21, D22) and the
## per-port delay tau (whole samples) from a discrete InternalDelay ss model,
## by reusing @ss/__ss_ext_build__ (the same helper c2d/d2c use) purely for
## data extraction.
function [B2, C2, D12, D21, D22, tau] = __internaldelay_ports__ (sys_dt)
  [ext, nu, ny] = __ss_ext_build__ (sys_dt);
  [~, Be, Ce, De] = ssdata (ext);
  B2  = Be(:, nu+1:end);
  C2  = Ce(ny+1:end, :);
  D12 = De(1:ny, nu+1:end);
  D21 = De(ny+1:end, 1:nu);
  D22 = De(ny+1:end, nu+1:end);
  tau = get (sys_dt, "internaldelay")(:);
endfunction


## Delay-buffer lookup: w(port) = z_hist(k - tau(port), port), read as 0 when
## k - tau(port) < 1 (no past sample yet -- the zero-history boundary).  At
## k == tau this still reads 0 (index 0); the first real-history read is at
## k == tau + 1, returning z_hist(1) -- i.e. exactly a tau-sample delay,
## matching the integer-sample shift used by __apply_timeresp_delay__.
##
## Assumes tau(port) >= 1 for every port: a port whose delay rounds to 0
## samples (an internal delay smaller than half a sample time) would read
## w(k) = z_hist(k), which is always 0 at read time since z_hist(k) is only
## written later in the same step -- silently dropping that port's D12/D22
## feedthrough instead of solving the resulting algebraic w=z loop. Not
## reachable via a well-formed nonzero delay; c2d's rounding is expected to
## keep tau >= 1 for any InternalDelay this function is asked to simulate.
function w = __delay_lookup__ (z_hist, k, tau, nports)
  w = zeros (nports, 1);
  for pp = 1 : nports
    idx = k - tau(pp);
    if (idx >= 1)
      w(pp) = z_hist(idx, pp);
    endif
  endfor
endfunction


function [tfinal, dt] = __sim_horizon__ (sys, tfinal, Ts)

  ## code based on __stepimp__.m of Kai P. Mueller and A. Scottedward Hodel

  N_MIN = 100;                                          # min number of points
  N_MAX = 10000;                                        # max number of points
  N_DEF = 2000;                                         # default number of points
  T_DEF = 10;                                           # default simulation time

  ev = pole (sys);

  TOL = max (abs (ev))*1.0e-10 + 2*eps;                 # values below TOL are assumed to be zero,
                                                        # avoid TOL = 0
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
  TOL = max (abs (ev))*1.0e-10 + 2*eps;                 # new TOL after bilinear transformation
  endif

  ## remove poles near zero from eigenvalue array ev
  nk = n;
    for k = 1 : n
    if (abs (ev(k)) < TOL)
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

    ev_conj_compl = ev (find (imag (ev) > TOL));
    if length (ev_conj_compl) > 0
      Dw0 = -real (ev_conj_compl);
      w   = imag (ev_conj_compl);
      t_osc = min ( 12*pi./w, 4./Dw0 );
      t_max_osc = max (t_osc);  # get max display time for slowest oscillation
    else
      t_max_osc = 0;
    endif

    if (continuous)
      dt = 0.1 * pi / ev_max;
    endif

    auto_tfinal = 0;  % flag for computed or given tfinal

    if (isempty (tfinal))
      ev_min = min (abs (ev));
      ev_real_min = min (abs (real (ev)));

      den = min ([ev_min, ev_real_min]);
      if (den < TOL)
        den =  max([ev_min, ev_real_min]);
      endif

      tfinal = 10 / den;
      auto_tfinal = 1;  # remeber that tfinal was computed, not given by the user

      tfinal = max (tfinal, t_max_osc); # make sure to show enough oscillations

      ## round up
      yy = 10^(ceil (log10 (tfinal)) - 1);
      tfinal = yy * ceil (tfinal / yy);
    endif

    if (continuous)

      ## Always select N such that tfinal < N*dt =< tfinal+dt
      N = fix (tfinal / dt) + 1;

      ## Ensure that tfinal is an integer multiple of dt and by
      ## the selection of N as above, we alwys reduce dt a little bit
      dt = tfinal/N;

      if (N < N_MIN)
        dt = tfinal / N_MIN;
      endif

      if (N > N_MAX)
        ## N is larger then N_MAX -> increase dt or reduce tfinal
        if (auto_tfinal)
          ## tfinal was computed: make it shorter and leave dt as it is in order to
          ## avoid aliasing
          tfinal = dt * N_MAX;  # adapt tfinal, not dt
          yy = 10^(ceil (log10 (tfinal)) - 1);  # round up again, since tfinal has changed
          tfinal = yy * ceil (tfinal / yy);
        else
          ## tfinal was selected by the user, do not change it, increase dt instead
          dt = tfinal / N_MAX;
        endif
      endif

    endif

  endif

  if (continuous && ! isempty (Ts))                     # catch case cont. system with dt specified
    dt = Ts;
  endif

endfunction

%!test  # no InternalDelay: unaffected (regression)
%! sys = ss (-1, 1, 1, 0);
%! y = step (sys);
%! assert (isempty (y), false);

%!test  # InternalDelay step vs an independent hand-written delay-buffer reference
%! T = 0.3; dt = 0.1;
%! sysd = c2d (ss (feedback (ss (-1, 1, 1, 0, "IODelay", T))), dt);
%! t = (0:dt:5)';
%! [y, ty] = step (sysd, t);
%! ## Independent reference: extract matrices (data plumbing only), simulate by hand.
%! [A, B1, C1, D11] = ssdata (sysd);
%! [ext, nu, ny] = __ss_ext_build__ (sysd);
%! [Ae, Be, Ce, De] = ssdata (ext);
%! B2 = Be(:, nu+1:end); C2 = Ce(ny+1:end, :);
%! D12 = De(1:ny, nu+1:end); D21 = De(ny+1:end, 1:nu); D22 = De(ny+1:end, nu+1:end);
%! tau = get (sysd, "internaldelay");
%! N = numel (ty); nst = rows (A);
%! xr = zeros (nst, 1); zh = zeros (N, 1); yref = zeros (N, 1);
%! for k = 1:N
%!   if (k - tau >= 1), w = zh(k - tau); else, w = 0; endif
%!   yref(k) = C1*xr + D11*1 + D12*w;
%!   zh(k) = C2*xr + D21*1 + D22*w;
%!   xr = A*xr + B1*1 + B2*w;
%! endfor
%! assert (y(:), yref, 1e-10);

%!test  # InternalDelay initial with nonzero x0: independent free-response reference
%! T = 0.3; dt = 0.1;
%! sysd = c2d (ss (feedback (ss (-1, 1, 1, 0, "IODelay", T))), dt);
%! x0 = 0.9;
%! t = (0:dt:5)';
%! [y, ty] = initial (sysd, x0, t);
%! [A, B1, C1, D11] = ssdata (sysd);
%! [ext, nu, ny] = __ss_ext_build__ (sysd);
%! [Ae, Be, Ce, De] = ssdata (ext);
%! B2 = Be(:, nu+1:end); C2 = Ce(ny+1:end, :);
%! D12 = De(1:ny, nu+1:end); D22 = De(ny+1:end, nu+1:end);
%! tau = get (sysd, "internaldelay");
%! N = numel (ty); nst = rows (A);
%! xr = x0(:); zh = zeros (N, 1); yref = zeros (N, 1);
%! for k = 1:N
%!   if (k - tau >= 1), w = zh(k - tau); else, w = 0; endif
%!   yref(k) = C1*xr + D12*w;
%!   zh(k) = C2*xr + D22*w;
%!   xr = A*xr + B2*w;
%! endfor
%! assert (y(:), yref, 1e-10);

%!test  # no InternalDelay: unaffected (regression)
%! sys = ss (-1, 1, 1, 0);
%! assert (isempty (impulse (sys)), false);
