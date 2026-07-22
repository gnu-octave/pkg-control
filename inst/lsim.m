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
## @deftypefn{Function File} {} lsim (@var{sys}, @var{u})
## @deftypefnx{Function File} {} lsim (@var{sys1}, @var{sys2}, ..., @var{sysN}, @var{u})
## @deftypefnx{Function File} {} lsim (@var{sys1}, @var{style1}, ..., @var{sysN}, @var{styleN}, @var{u})
## @deftypefnx{Function File} {} lsim (@var{sys1}, ..., @var{u}, @var{t})
## @deftypefnx{Function File} {} lsim (@var{sys1}, ..., @var{u}, @var{t}, @var{x0})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} lsim (@var{sys}, @var{u})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} lsim (@var{sys}, @var{u}, @var{t})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} lsim (@var{sys}, @var{u}, @var{t}, @var{x0})
## @deftypefnx{Function File} {[...] =} lsim (..., @var{method})
##
## Simulate @acronym{LTI} model response to arbitrary inputs.  If no output arguments are given,
## the system response is plotted on the screen.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} model.  System must be proper, i.e. it must not have more zeros than poles.
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
## Vector of initial conditions for each state of a system in state space.
## If not specified, a zero vector is assumed.
## Note: A vector @var{x0} provided for an input-output system representation
## is ignored and a zero vector of initial conditions is used instead because
## the internally used state space representaton does generally not match
## the one assumed for @var{x0}. For a simulation of an input-output model
## with initial conditions for the output @math{y} and its time-derivatives,
## see remarks below.
## @item style
## Line style and color, e.g. 'r' for a solid red line or '-.k' for a dash-dotted
## black line. See @command{help plot} for details.
## @item method
## Method that is used to discretize a continuous-time system for the simulation.
## See @ref{@@lti/c2d} for possible methods. If @var{method} is not provided,
## the default is 'foh'.
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
## @strong{Remarks}
##
## @itemize
##
## @item For the simulation, continuous-time systems are discretized using the
## selected method @var{method} or the default first-order-hold method (foh), see
## @ref{@@lti/c2d} for details.
##
## @item For a SISO input-output model @math{G} and initial values for the output
## @math{y} and its derivatives up to order @math{n-1} the corresponding state space
## represetaiton is computed by:
##
## @example
## @group
## [A,b,c,d] = ssdata (G);
## T = obsv (A, c);
## G_ss = ss2ss (ss (G), T);
## initial (G_ss, x0);
## @end group
## @end example
##
## Note that, in general, the states of @math{G_ss} are only equal to the output
## @math{y} and its first @math{n-1} time derivaties if @math{u=0}, which is the
## case for the initial conditions immediately before @math{t=0}.
##
## @end itemize
##
## @seealso{impulse, initial, step}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.5

function [y_r, t_r, x_r] = lsim (varargin)

  ## TODO: individual initial state vectors 'x0' for each system
  ##       there would be conflicts with other arguments,
  ##       maybe a cell {x0_1, x0_2, ..., x0_N} would be a solution?

  if (nargin < 2)
    print_usage ();
  endif

  ## if method is provided, it is the last parameter and method
  ## is the only char parameter at the end of the parameter list
  if (ischar (varargin{end}))
    method = varargin{end};
    varargin = varargin(1:end-1);
  else
    method = "foh";
  endif

  ## get remaining parameters

  ## FIXME: what is the prupose of the following code, which converts
  ##        logical args into double. Where an logical args come from?
  ## idx = cellfun (@islogical, varargin);
  ## tmp = cellfun (@double, varargin(idx), "uniformoutput", false);
  ## varargin(idx) = tmp;

  names = arrayfun (@inputname, 1:nargin, 'uniformoutput', false);

  [sys_idx, sty_idx, mat_idx, names, sty] = __control_args__ (varargin, {"@lti", @ischar, @is_real_matrix}, names);

  inv_idx = ! (sys_idx | mat_idx | sty_idx);            # invalid arguments

  if (any (inv_idx))
    warning ("lsim: arguments number %s are invalid and are being ignored\n", ...
             mat2str (find (inv_idx)(:).'));
  endif

  if (nnz (sys_idx) == 0)
    error ("lsim: require at least one LTI model");
  endif

  if (nargout > 0 && (nnz (sys_idx) > 1 || any (sty_idx)))
    print_usage ();
  endif

  if (! size_equal (varargin{sys_idx}))
    error ("lsim: all LTI models must have equal size");
  endif

  if (any (find (sty_idx) < find (sys_idx)(1)))
    warning ("lsim: strings in front of first LTI model are being ignored\n");
  endif

  t = [];  x0 = [];                                     # default arguments

  switch (nnz (mat_idx))
    case 0
      error ("lsim: require input signal 'u'");
    case 1
      u = varargin{mat_idx};
    case 2
      [u, t] = varargin{mat_idx};
    case 3
      [u, t, x0] = varargin{mat_idx};
    otherwise
      print_usage ();
  endswitch

  if (is_real_vector (u))                               # allow row vectors for single-input systems
    u = vec (u);
  elseif (isempty (u))                                  # ! is_real_matrix (u)  already tested
    error ("lsim: input signal 'u' must be a real-valued matrix");
  endif

  if (! is_real_vector (t) && ! isempty (t))
    error ("lsim: time vector 't' must be real-valued or empty");
  endif

  if (! isequal (t, unique (t)))
    error ("lsim: time vector 't' must be sorted");
  endif

  if (! is_real_vector (x0) && ! isempty (x0))
    error ("lsim: initial state vector 'x0' must be empty or a real-valued vector");
  endif

  n_sys = nnz (sys_idx);

  ## function [y, t, x_arr] = __linear_simulation__ (sys, u, t, x0)

  if (! isempty (x0))
    is_ss = cellfun (@isa, varargin(sys_idx), {"ss"});
    if (! all (is_ss))
      no_ss = find (is_ss==0);
      no_ss_names = sprintf ("%s, ", names{no_ss});
      no_ss_names = no_ss_names(1:end-2);
      warning ("lsim: system %s not in state space, x0 is ambiguous and is ignored assuming zero initial conditions\n", no_ss_names);
      for j = 1:length (no_ss)
        names{j} = [names{j}, " (x0=0)"];
      endfor
    endif
  endif

  [y, t, x] = cellfun (@__linear_simulation__, varargin(sys_idx), {u}, {t}, {x0}, {method}, "uniformoutput", false);


  if (nargout == 0)                                     # plot information

    [p, m] = size (varargin(sys_idx){1});
    ct_idx = cellfun (@isct, varargin(sys_idx));
    str = "Linear Simulation Results";

    outname = get (varargin(sys_idx){end}, "outname");
    outname = __labels__ (outname, "y");

    inname = get (varargin(sys_idx){end}, "inname");
    inname = __labels__ (inname, "u");
    incolororder = 0.6 * ones (m,3);
    for j = 2:m
      incolororder(j,:) = incolororder(j-1,:) + (j-1)*0.2/(m-1);
    endfor

    for k = 1 : n_sys                                   # for every system
      if (ct_idx(k))                                    # continuous-time system
        for i = 1 : p                                   # for every output
          if (p != 1)
            subplot (p, 1, i);
          endif
          plot (t{k}, y{k}(:, i), sty{k}{:});
          hold on;
          grid on;
          if (k == n_sys)
            ylabel (outname{i});
            if (i == 1)
              title (str);
            endif
            set (gca (), "colororder", incolororder);
            plot (t{k}, u, '-.');  # plot input
            axis tight
            ylim (__axis_margin__ (ylim))
            legend ([names,inname']);
          endif
        endfor
      else                                              # discrete-time system
        for i = 1 : p                                   # for every output
          if (p != 1)
            subplot (p, 1, i);
          endif
          stairs (t{k}, y{k}(:, i), sty{k}{:});
          hold on;
          grid on;
          if (k == n_sys)
            ylabel (outname{i});
            if (i == 1)
              title (str);
            endif
            set (gca (), "colororder", incolororder);
            plot (t{k}, u, '-.');  # plot input
            axis tight
            ylim (__axis_margin__ (ylim))
            legend ([names,inname']);
          endif
        endfor
      endif
    endfor
    xlabel ("Time [s]");
%    if (p == 1 && m == 1)
%      u1_name = __labels__ (get (varargin(sys_idx){1}, "inname"), "u");
%      names{end+1} = u1_name{1};
%      legend (names)
%    endif
    hold off;
  else                                                  # return values
    y_r = y{1};
    t_r = t{1};
    x_r = x{1};
  endif

endfunction


function [y, t, x_arr] = __linear_simulation__ (sys, u, t, x0, method)

  if (! isa (sys, "ss"))
    x0 =[]; # ignore initial condition for system not in state space
  endif

  [urows, ucols] = size (u);
  len_t = length (t);

  if (isct (sys))                               # continuous-time system
    was_ct = 1;
    if (isempty (t))                            # lsim (sys, u, [], ...)
      error ("lsim: time vector 't' must not be empty");
    elseif (len_t == 1)                         # lsim (sys, u, tfinal, ...)
      dt = t / (urows - 1);
      t = vec (linspace (0, t, urows));
    elseif ((len_t != urows) && (len_t != ucols))
      error ("lsim: length of time vector (%d) doesn't match input signal (%dx%d) or (%dx%d)\n", ...
             len_t, urows, ucols, ucols, urows);
    else                                        # lsim (sys, u, t, ...)
      if (len_t == ucols)
        u = u';
        [urows, ucols] = size (u);
      endif
      dt = abs (t(end) - t(1)) / (urows - 1);   # assume that t is regularly spaced
      t = vec (linspace (t(1), t(end), urows));
    endif
    sys = c2d (ss (sys), dt, method);           # convert to discrete-time model (in ss for accuracy)
  else                                          # discrete-time system
    was_ct = 0;
    dt = abs (get (sys, "tsam"));               # use 1 second as default if tsam is unspecified (-1)
    if (isempty (t))                            # lsim (sys, u)
      m = length (sys.inputname);               # we can not verify shape of u by length t
      if ((ucols != m) && (urows != m))
        error ("lsim: input vector 'u' must have %d columns or rows", m);
      else
        if (urows == m)
          u = u';
          [urows, ucols] = size (u);
        endif
      endif
      t = vec (linspace (0, dt*(urows-1), urows));
    elseif (len_t == 1)                         # lsim (sys, u, tfinal)
      ## TODO: maybe raise warning if  abs (tfinal - dt*(urows-1)) > dt
      t = vec (linspace (0, dt*(urows-1), urows));
    elseif ((len_t != urows) && (len_t != ucols))
      error ("lsim: length of time vector (%d) doesn't match input signal (%dx%d) or (%dx%d)\n", ...
             len_t, urows, ucols, ucols, urows);
      if (len_t == ucols)
        u = u';
        [urows, ucols] = size (u);
      endif
    else                                        # lsim (sys, u, t, ...)
      t = vec (linspace (t(1), t(end), len_t));
    endif
  endif

  [A, B, C, D] = ssdata (sys);
  [p, m] = size (D);                            # number of outputs and inputs
  n = rows (A);                                 # number of states

  if (ucols != m)
    error ("lsim: input vector 'u' must have %d columns", m);
  endif

  if (isempty (x0))
    x0 = zeros (n, 1);
  elseif (n != length (x0) || ! is_real_vector (x0))
    error ("lsim: 'x0' must be a vector with %d elements", n);
  endif

  x0 = vec (x0);                                # make sure that x0 is a column vector

  is_foh = was_ct && strcmp (method, "foh") && (max (size (sys.userdata)) > 0);

  if (hasinternaldelay (sys))

    ## InternalDelay: the internal-delay port input w(k) is drawn from a
    ## buffer of already-computed past z values (tau samples ago), so there
    ## is no algebraic loop -- just a per-step buffer lookup.  The recursion
    ## is linear in (x0, u), so it composes with x0 by direct state basis (no
    ## superposition trick needed) and, when ordinary I/O delay is ALSO
    ## present, with the prior phase's post-hoc __apply_timeresp_delay__ shift
    ## via the same per-channel superposition it already uses.  (An
    ## InternalDelay system is discretized through Task 1's c2d path, which
    ## never populates sys.userdata, so is_foh is always false here.)
    [B2, C2, D12, D21, D22, tau] = __internaldelay_ports__ (sys);

    if (! hasdelay (sys))
      ## pure internal delay: one combined buffered recursion (x0 and u together)
      [y, x_arr] = __buffered_sim__ (A, B, C, D, B2, C2, D12, D21, D22, tau, ...
                                     x0, u, urows, p, n);
    else
      ## internal delay AND ordinary I/O delay: superpose a zero-input term
      ## (x0 alone, shifted by OutputDelay) with one zero-state term per input
      ## channel (shifted by that channel's totaldelay); each sub-simulation
      ## carries its own internal-delay buffer, and summing reconstructs the
      ## full response because the buffer recursion is linear.
      total = totaldelay (sys);            # p-by-m, whole samples
      outdelay = get (sys, "outputdelay"); # p-by-1, whole samples

      y = zeros (urows, p);
      x_arr = zeros (urows, n);

      [y_free, x_free] = __buffered_sim__ (A, B, C, D, B2, C2, D12, D21, D22, ...
                                           tau, x0, zeros (urows, m), urows, p, n);
      x_arr += x_free;
      for row = 1 : p
        y(:, row) += __apply_timeresp_delay__ (y_free(:, row), outdelay(row));
      endfor

      for j = 1 : m
        uj = zeros (urows, m);
        uj(:, j) = u(:, j);
        [y_j, x_j] = __buffered_sim__ (A, B, C, D, B2, C2, D12, D21, D22, ...
                                       tau, zeros (n, 1), uj, urows, p, n);
        x_arr += x_j;
        for row = 1 : p
          y(:, row) += __apply_timeresp_delay__ (y_j(:, row), total(row, j));
        endfor
      endfor
    endif

  elseif (! hasdelay (sys))

    ## preallocate memory
    y = zeros (urows, p);
    x_arr = zeros (urows, n);

    x = x0;

    ## When discretization method was foh transform initial state into
    ## the states representing the foh form.
    ## The required matrix "Bd1" is stored by c2d in sys.userdata
    if (is_foh)
      x = x - sys.userdata * u(1,:)';
    endif

    ## simulation
    for k = 1 : urows
      y(k, :) = C * x  +  D * u(k, :).';
      x_arr(k, :) = x;
      x = A * x  +  B * u(k, :).';
    endfor

    ## When discretization method was foh transform back from foh states
    ## into original state
    if (is_foh)
      x_arr = x_arr + u * sys.userdata';
    endif

  else

    ## Delayed system: superposition of a zero-input term (from x0 alone,
    ## only OutputDelay applies -- there is no input for InputDelay/IODelay
    ## to act on) plus one zero-state term per input channel (from u(:,j)
    ## alone, shifted by that channel's full totaldelay, which already
    ## includes OutputDelay). Both terms use the same delay-free A,B,C,D,
    ## so their states add linearly and are left unshifted (see plan/spec).
    total = totaldelay (sys);            # p-by-m, whole samples
    outdelay = get (sys, "outputdelay"); # p-by-1, whole samples

    y = zeros (urows, p);
    x_arr = zeros (urows, n);

    ## zero-input term: driven by x0 alone, with no associated input
    ## sample, so no foh initial-state correction applies here -- the
    ## foh correction for u(1,:) is fully accounted for by the per-column
    ## offsets below (each column contributes its own -Bd1(:,j)*u(1,j)
    ## piece); applying it here too would double-count it.
    x = x0;
    y_free = zeros (urows, p);
    for k = 1 : urows
      y_free(k, :) = C * x;
      x_arr(k, :) = x;
      x = A * x;
    endfor
    for row = 1 : p
      y(:, row) += __apply_timeresp_delay__ (y_free(:, row), outdelay(row));
    endfor

    ## zero-state term, one input channel at a time
    for j = 1 : m
      uj = zeros (urows, m);
      uj(:, j) = u(:, j);

      xj = zeros (n, 1);
      if (is_foh)
        xj = -sys.userdata(:, j) * u(1, j);
      endif

      y_j = zeros (urows, p);
      x_arr_j = zeros (urows, n);
      for k = 1 : urows
        y_j(k, :) = C * xj  +  D * uj(k, :).';
        x_arr_j(k, :) = xj;
        xj = A * xj  +  B * uj(k, :).';
      endfor
      if (is_foh)
        x_arr_j = x_arr_j + uj * sys.userdata';
      endif

      x_arr += x_arr_j;
      for row = 1 : p
        y(:, row) += __apply_timeresp_delay__ (y_j(:, row), total(row, j));
      endfor
    endfor

  endif

endfunction


## Extract the internal-delay port matrices (B2, C2, D12, D21, D22) and the
## per-port delay tau (whole samples) from a discrete InternalDelay ss model,
## reusing @ss/__ss_ext_build__ (the same helper c2d/d2c use) for data
## extraction only.
function [B2, C2, D12, D21, D22, tau] = __internaldelay_ports__ (sys)
  [ext, nu, ny] = __ss_ext_build__ (sys);
  [~, Be, Ce, De] = ssdata (ext);
  B2  = Be(:, nu+1:end);
  C2  = Ce(ny+1:end, :);
  D12 = De(1:ny, nu+1:end);
  D21 = De(ny+1:end, 1:nu);
  D22 = De(ny+1:end, nu+1:end);
  tau = get (sys, "internaldelay")(:);
endfunction


## One buffered per-sample simulation of an InternalDelay ss model over the
## whole horizon: at each step the internal-delay port input w is read from a
## buffer of past z values (tau samples ago; zero when k - tau < 1 -- the
## zero-history boundary), never the current step's z.  Returns the output y
## and state trajectory x_arr; linear in (x_init, umat).
##
## Assumes tau(port) >= 1 for every port, same as __delay_lookup__ in
## __time_response__.m: a port whose delay rounds to 0 samples would read
## w(k) = z_hist(k) before it is written this step, silently dropping that
## port's D12/D22 feedthrough instead of solving the resulting algebraic
## w=z loop.  Not reachable via a well-formed nonzero InternalDelay.
function [y, x_arr] = __buffered_sim__ (A, B, C, D, B2, C2, D12, D21, D22, ...
                                        tau, x_init, umat, urows, p, n)
  nports = numel (tau);
  z_hist = zeros (urows, nports);
  y = zeros (urows, p);
  x_arr = zeros (urows, n);
  x = x_init;
  for k = 1 : urows
    w = zeros (nports, 1);
    for pp = 1 : nports
      idx = k - tau(pp);
      if (idx >= 1)
        w(pp) = z_hist(idx, pp);
      endif
    endfor
    uk = umat(k, :).';
    y(k, :) = C * x + D * uk + D12 * w;
    z_hist(k, :) = (C2 * x + D21 * uk + D22 * w).';
    x_arr(k, :) = x;
    x = A * x + B * uk + B2 * w;
  endfor
endfunction


%!test
%! n = 5;
%! m = 3;
%! p = 2;
%! A = diag ([0:-2:-2*(n-1)]);
%! B = [ (1:1:n)' (-1:1:n-2)' (2:1:n+1)'];
%! C = [1 0 1 0 0 ; 0 1 0 1 1 ];
%! D = zeros (p,m);
%!
%! sys = ss(A, B, C, D);
%! dt = 0.1;
%! t = 0:dt:1;
%! x0 = zeros(n,1);
%! u = [ sin(2*t') cos(3*t') sin(4*t') ];
%! [y1, t1] = lsim(sys, u, t, x0);
%! [y2, t2] = lsim(sys, u', t, x0);
%!
%! sysd = c2d (sys, dt, 'foh');
%! x0 = x0 - sysd.userdata * u(1,:)'; # foh-doscretization
%! [y3, t3] = lsim(sysd, u, [], x0);
%! [y4, t4] = lsim(sysd, u', [], x0);
%!
%! assert (y1,y2,1e-4);
%! assert (y1,y3,1e-4);
%! assert (y1,y4,1e-4);

%!test  # SISO InputDelay: lsim output matches a manually-shifted delay-free lsim
%! sys = tf (1, [1 1], "InputDelay", 0.3);
%! sys_nodelay = tf (1, [1 1]);
%! t = 0:0.05:5;
%! u = sin (t)';
%! [y, tt] = lsim (sys, u, t);
%! [y0, tt0] = lsim (sys_nodelay, u, t);
%! dt = tt(2) - tt(1);
%! k = round (0.3 / dt);
%! expected = [zeros(k, 1); y0(1:end-k)];
%! assert (y, expected, 1e-6);

%!test  # MIMO IODelay: each output channel is the sum of per-input contributions,
%!       # each shifted by its own total delay
%! sys = tf ({1, 1; 1, 1}, {[1 1], [1 2]; [1 3], [1 4]});
%! sys = set (sys, "IODelay", [0.2, 0; 0, 0.4]);
%! sys_nodelay = tf ({1, 1; 1, 1}, {[1 1], [1 2]; [1 3], [1 4]});
%! t = (0:0.05:5)';
%! u = [sin(t), cos(t)];
%! [y, tt] = lsim (sys, u, t);
%! [y0, tt0] = lsim (sys_nodelay, u, t);
%! dt = tt(2) - tt(1);
%! total = [0.2, 0; 0, 0.4];
%! ## reconstruct expected per-output response by re-simulating each input alone
%! u1 = [u(:,1), zeros(size(t))];
%! u2 = [zeros(size(t)), u(:,2)];
%! [y01] = lsim (sys_nodelay, u1, t);
%! [y02] = lsim (sys_nodelay, u2, t);
%! k11 = round (total(1,1)/dt); k12 = round (total(1,2)/dt);
%! k21 = round (total(2,1)/dt); k22 = round (total(2,2)/dt);
%! if (k11 == 0)
%!   s11 = y01(:,1);
%! else
%!   s11 = [zeros(k11,1); y01(1:end-k11,1)];
%! endif
%! if (k12 == 0)
%!   s12 = y02(:,1);
%! else
%!   s12 = [zeros(k12,1); y02(1:end-k12,1)];
%! endif
%! if (k21 == 0)
%!   s21 = y01(:,2);
%! else
%!   s21 = [zeros(k21,1); y01(1:end-k21,2)];
%! endif
%! if (k22 == 0)
%!   s22 = y02(:,2);
%! else
%!   s22 = [zeros(k22,1); y02(1:end-k22,2)];
%! endif
%! expected1 = s11 + s12;
%! expected2 = s21 + s22;
%! assert (y(:,1), expected1, 1e-6);
%! assert (y(:,2), expected2, 1e-6);

%!test  # nonzero x0 on a delayed ss system: zero-input term shifted by OutputDelay only
%! A = -1; B = 1; C = 1; D = 0;
%! sys = ss (A, B, C, D, "OutputDelay", 0.3);
%! sys_nodelay = ss (A, B, C, D);
%! t = (0:0.05:5)';
%! u = zeros (size (t));       # isolate the zero-input term
%! x0 = 2;
%! [y, tt] = lsim (sys, u, t, x0);
%! [y0, tt0] = lsim (sys_nodelay, u, t, x0);
%! dt = tt(2) - tt(1);
%! k = round (0.3 / dt);
%! expected = [zeros(k, 1); y0(1:end-k)];
%! assert (y, expected, 1e-6);

%!demo
%! clf;
%! A = [-3   0   0;
%!       0  -2   1;
%!      10 -17   0];
%! B = [4  0;
%!      0 -1;
%!      0 -1];
%! C = [0 0 1;
%!      1 2 0];
%! D = [ 0 0;
%!       0 0 ];
%! S = ss(A,B,C,D);
%! t = 0:0.01:6;
%! u = [ 0.2+0.3*sin(1.3*t') , cos(2*t') ];
%! x0 = [0 0.1 0];
%! lsim(S, u, t, x0);

%!test  # InternalDelay lsim vs an independent hand-written delay-buffer reference
%! T = 0.3; dt = 0.1;
%! sysd = c2d (ss (feedback (ss (-1, 1, 1, 0, "IODelay", T))), dt);
%! t = (0:dt:3)';
%! u = sin (t);
%! [y, tt, x] = lsim (sysd, u, t);
%! ## Independent reference: extract the matrices (data plumbing only) and run
%! ## the delay-buffer recursion by hand with plain array indexing.
%! [A, B1, C1, D11] = ssdata (sysd);
%! [ext, nu, ny] = __ss_ext_build__ (sysd);
%! [Ae, Be, Ce, De] = ssdata (ext);
%! B2 = Be(:, nu+1:end); C2 = Ce(ny+1:end, :);
%! D12 = De(1:ny, nu+1:end); D21 = De(ny+1:end, 1:nu); D22 = De(ny+1:end, nu+1:end);
%! tau = get (sysd, "internaldelay");
%! N = numel (t); nst = rows (A);
%! xr = zeros (nst, 1); zh = zeros (N, 1); yref = zeros (N, 1); xref = zeros (N, nst);
%! for k = 1:N
%!   if (k - tau >= 1), w = zh(k - tau); else, w = 0; endif
%!   yref(k) = C1*xr + D11*u(k) + D12*w;
%!   zh(k) = C2*xr + D21*u(k) + D22*w;
%!   xref(k, :) = xr.';
%!   xr = A*xr + B1*u(k) + B2*w;
%! endfor
%! assert (y, yref, 1e-10);
%! assert (x, xref, 1e-10);

%!test  # nonzero x0 composes by plain linear superposition (direct state basis)
%! T = 0.3; dt = 0.1;
%! sysd = c2d (ss (feedback (ss (-1, 1, 1, 0, "IODelay", T))), dt);
%! t = (0:dt:3)'; u = sin (t); x0 = 0.7;
%! y_full = lsim (sysd, u, t, x0);
%! y_zero = lsim (sysd, u, t, 0);
%! y_free = lsim (sysd, zeros (size (t)), t, x0);
%! assert (y_full, y_zero + y_free, 1e-10);
