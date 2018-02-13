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
## @deftypefnx{Function File} {} lsim (@var{sys1}, @var{sys2}, @dots{}, @var{sysN}, @var{u})
## @deftypefnx{Function File} {} lsim (@var{sys1}, @var{'style1'}, @dots{}, @var{sysN}, @var{'styleN'}, @var{u})
## @deftypefnx{Function File} {} lsim (@var{sys1}, @dots{}, @var{u}, @var{t})
## @deftypefnx{Function File} {} lsim (@var{sys1}, @dots{}, @var{u}, @var{t}, @var{x0})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} lsim (@var{sys}, @var{u})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} lsim (@var{sys}, @var{u}, @var{t})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} lsim (@var{sys}, @var{u}, @var{t}, @var{x0})
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
## Vector of initial conditions for each state.  If not specified, a zero vector is assumed.
## @item 'style'
## Line style and color, e.g. 'r' for a solid red line or '-.k' for a dash-dotted
## black line.  See @command{help plot} for details.
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
## Version: 0.5

function [y_r, t_r, x_r] = lsim (varargin)

  ## TODO: individual initial state vectors 'x0' for each system
  ##       there would be conflicts with other arguments,
  ##       maybe a cell {x0_1, x0_2, ..., x0_N} would be a solution?

  if (nargin < 2)
    print_usage ();
  endif

  idx = cellfun (@islogical, varargin);
  tmp = cellfun (@double, varargin(idx), "uniformoutput", false);
  varargin(idx) = tmp;

  sys_idx = cellfun (@isa, varargin, {"lti"});          # LTI models
  mat_idx = cellfun (@is_real_matrix, varargin);        # matrices
  sty_idx = cellfun (@ischar, varargin);                # string (style arguments)

  inv_idx = ! (sys_idx | mat_idx | sty_idx);            # invalid arguments

  if (any (inv_idx))
    warning ("lsim: arguments number %s are invalid and are being ignored", ...
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
    warning ("lsim: strings in front of first LTI model are being ignored");
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


  ## function [y, t, x_arr] = __linear_simulation__ (sys, u, t, x0)
  
  [y, t, x] = cellfun (@__linear_simulation__, varargin(sys_idx), {u}, {t}, {x0}, "uniformoutput", false);


  if (nargout == 0)                                     # plot information
    ## extract plotting styles
    tmp = cumsum (sys_idx);
    tmp(sys_idx | ! sty_idx) = 0;
    n_sys = nnz (sys_idx);
    sty = arrayfun (@(x) varargin(tmp == x), 1:n_sys, "uniformoutput", false);

    ## default plotting styles if empty
    colororder = get (gca, "colororder");
    rc = rows (colororder);
    def = arrayfun (@(k) {"color", colororder(1+rem (k-1, rc), :)}, 1:n_sys, "uniformoutput", false);
    idx = cellfun (@isempty, sty);
    sty(idx) = def(idx);
  
    ## get system names for legend
    ## leg = cellfun (@inputname, find (sys_idx), "uniformoutput", false);
    leg = cell (1, n_sys);
    idx = find (sys_idx);
    for k = 1 : n_sys
      try
        leg{k} = inputname (idx(k));
      catch
        leg{k} = "";                                    # catch case  lsim (lticell{:}, ...)
      end_try_catch
    endfor
    
    [p, m] = size (varargin(sys_idx){1});
    ct_idx = cellfun (@isct, varargin(sys_idx));
    str = "Linear Simulation Results";
    outname = get (varargin(sys_idx){end}, "outname");
    outname = __labels__ (outname, "y");

    for k = 1 : n_sys                                   # for every system
      if (ct_idx(k))                                    # continuous-time system                                           
        for i = 1 : p                                   # for every output
          if (p != 1)
            subplot (p, 1, i);
          endif
          plot (t{k}, y{k}(:, i), sty{k}{:});
          hold on;
          # input should be plotted in the background using uistack, which isn't
          # implemented yet
          plot (t{k}, u, 'Color', [0.5 0.5 0.5]);       # plot input
          grid on;
          if (k == n_sys)
            axis tight
            ylim (__axis_margin__ (ylim))
            ylabel (outname{i});
            if (i == 1)
              title (str);
            endif
          endif
        endfor
      else                                              # discrete-time system
        for i = 1 : p                                   # for every output
          if (p != 1)
            subplot (p, 1, i);
          endif
          stairs (t{k}, y{k}(:, i), sty{k}{:});
          hold on;
          # input should be plotted in the background using uistack, which isn't
          # implemented yet
          plot (t{k}, u, 'Color', [0.5 0.5 0.5]);       # plot input
          grid on;
          if (k == n_sys)
            axis tight;
            ylim (__axis_margin__ (ylim))
            ylabel (outname{i});
            if (i == 1)
              title (str);
            endif
          endif
        endfor
      endif
    endfor
    xlabel ("Time [s]");
    if (p == 1 && m == 1)
      legend (leg)
    endif
    hold off;
  else                                                  # return values
    y_r = y{1};
    t_r = t{1};
    x_r = x{1};
  endif
  
endfunction


function [y, t, x_arr] = __linear_simulation__ (sys, u, t, x0)

  method = "zoh";
  [urows, ucols] = size (u);
  len_t = length (t);

  if (isct (sys))                               # continuous-time system
    if (isempty (t))                            # lsim (sys, u, [], ...)
      error ("lsim: time vector 't' must not be empty");
    elseif (len_t == 1)                         # lsim (sys, u, tfinal, ...)
      dt = t / (urows - 1);
      t = vec (linspace (0, t, urows));
    elseif (len_t != urows)
      error ("lsim: length of time vector (%d) doesn't match input signal (%dx%d)", ...
             len_t, urows, ucols);
    else                                        # lsim (sys, u, t, ...)
      dt = abs (t(end) - t(1)) / (urows - 1);   # assume that t is regularly spaced
      t = vec (linspace (t(1), t(end), urows));
    endif
    sys = c2d (sys, dt, method);                # convert to discrete-time model
  else                                          # discrete-time system
    dt = abs (get (sys, "tsam"));               # use 1 second as default if tsam is unspecified (-1)
    if (isempty (t))                            # lsim (sys, u)
      t = vec (linspace (0, dt*(urows-1), urows));
    elseif (len_t == 1)                         # lsim (sys, u, tfinal)
      ## TODO: maybe raise warning if  abs (tfinal - dt*(urows-1)) > dt
      t = vec (linspace (0, dt*(urows-1), urows));
    elseif (len_t != urows)
      error ("lsim: length of time vector (%d) doesn't match input signal (%dx%d)", ...
             len_t, urows, ucols);
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

  ## preallocate memory
  y = zeros (urows, p);
  x_arr = zeros (urows, n);

  ## initial conditions
  if (isempty (x0))
    x0 = zeros (n, 1);
  elseif (n != length (x0) || ! is_real_vector (x0))
    error ("lsim: 'x0' must be a vector with %d elements", n);
  endif

  x = vec (x0);                                 # make sure that x is a column vector

  ## simulation
  for k = 1 : urows
    y(k, :) = C * x  +  D * u(k, :).';
    x_arr(k, :) = x;
    x = A * x  +  B * u(k, :).';
  endfor

endfunction


## TODO: add test cases
