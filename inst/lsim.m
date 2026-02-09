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
## See
## @inlinefmtifelse{latex, @link{@@lti/c2d,@@lti/c2d}, @ref{@@lti/c2d}}
## for possible methods. If @var{method} is not provided, the default is 'foh'.
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
## @inlinefmtifelse{latex, @link{@@lti/c2d,@@lti/c2d}, @ref{@@lti/c2d}}
## for details.
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
  idx = cellfun (@islogical, varargin);
  tmp = cellfun (@double, varargin(idx), "uniformoutput", false);
  varargin(idx) = tmp;

  sys_idx = cellfun (@isa, varargin, {"lti"});          # LTI models
  mat_idx = cellfun (@is_real_matrix, varargin);        # matrices
  sty_idx = cellfun (@ischar, varargin);                # string (style arguments)

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

  ## get system names
  leg = cell (1, n_sys);
  idx = find (sys_idx);
  names = cell (1, n_sys);
  for k = 1 : n_sys
    try
      names{k} = inputname (idx(k));
    catch
      names{k} = ['Sys',num2str(k)];    # catch case  lsim (lticell{:}, ...)
    end_try_catch
  endfor

  ## function [y, t, x_arr] = __linear_simulation__ (sys, u, t, x0)

  if (! isempty (x0))
    is_ss = cellfun (@isa, varargin(sys_idx), {"ss"});
    if (! all (is_ss))
      no_ss = find (is_ss==0);
      no_ss_names = sprintf ("%s, ", names{no_ss});
      no_ss_names = no_ss_names(1:end-2);
      warning ("lsim: system %s not in state space, x0 is ambiguous and system is ignored\n", no_ss_names);
      for j = 1:length (no_ss)
        names{j} = [names{j}, " (x0=0)"];
      endfor
    endif
  endif

  [y, t, x] = cellfun (@__linear_simulation__, varargin(sys_idx), {u}, {t}, {x0}, {method}, "uniformoutput", false);


  if (nargout == 0)                                     # plot information
    ## extract plotting styles
    tmp = cumsum (sys_idx);
    tmp(sys_idx | ! sty_idx) = 0;
    sty = arrayfun (@(x) varargin(tmp == x), 1:n_sys, "uniformoutput", false);

    ## default plotting styles if empty
    colororder = get (gca, "colororder");
    rc = rows (colororder);
    def = arrayfun (@(k) {"color", colororder(1+rem (k-1, rc), :)}, 1:n_sys, "uniformoutput", false);
    idx = cellfun (@isempty, sty);
    sty(idx) = def(idx);

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
      legend (names)
    endif
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

  ## When discretization method was foh transform initial state into
  ## the states representing the foh form.
  ## The required matrix "Bd1" is stored by c2d in sys.userdata
  if was_ct && (strcmp (method, "foh")) && (max (size (sys.userdata)) > 0)
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
  if was_ct && (strcmp (method, "foh")) && (max (size (sys.userdata)) > 0)
    x_arr = x_arr + u * sys.userdata';
  endif

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

%!demo
%! clf;
%! A = [-3   0   0;
%!       0  -2   1;
%!      10 -17   0];
%! B = [4;
%!      0;
%!      0];
%! C = [0 0 1];
%! D = 0;
%! sys = ss(A,B,C,D);
%! t = 0:0.01:10;
%! u = zeros (length(t) ,1);
%! x0 = [0 0.1 0];
%! lsim(sys, u, t, x0);
