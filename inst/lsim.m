## Copyright (C) 2009   Lukas F. Reichlin
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
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn{Function File} {[@var{y}, @var{t}, @var{x}] =} lsim (@var{sys}, @var{u}, @var{t})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} lsim (@var{sys}, @var{u}, @var{t}, @var{x0})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} lsim (@var{sys}, @var{u}, @var{t}, @var{x0}, @var{method})
## Simulate LTI model response to arbitrary inputs. If no output arguments are given,
## the system response is plotted on the screen.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.1

function [y_r, t_r, x_r] = lsim (sys, u, t = [], x0 = [], method = "zoh")

  if (nargin < 2 || nargin > 5)
    print_usage ();
  endif

  if (! isa (sys, "ss"))
    sys = ss (sys);  # sys must be proper
  endif

  [A, B, C, D, tsam] = ssdata (sys);

  digital = ! isct (sys);  # static gains are treated as analog systems

  urows = rows (u);
  ucols = columns (u);

  if (digital)  # discrete system
    if (isempty (t))
      dt = tsam;
      tfinal = tsam * (urows - 1);
    elseif (length (t) == 1)
      dt = tsam;
      tfinal = t;
    else
      warning ("lsim: spacing of time vector has no effect on sampling time of digital system");
      dt = tsam;
      tfinal = t(end);
    endif
  else  # continuous system
    if (isempty (t))
      error ("lsim: invalid time vector");
    elseif (length (t) == 1)
      dt = t / (urows - 1);
      tfinal = t;
    else
      dt = t(2) - t(1);  # assume that t is regularly spaced
      tfinal = t(end);
    endif
    sys = c2d (sys, dt, method);  # convert to discrete model
  endif

  [F, G] = ssdata (sys);  # matrices C and D don't change

  n = rows (F);  # number of states
  m = columns (G);  # number of inputs
  p = rows (C);  # number of outputs

  t = (0 : dt : tfinal)';  # time vector
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
  elseif (n != length (x0))
    error ("initial: x0 must be a vector with %d elements", n);
  endif

  x = x0(:);  # make sure that x is a column vector

  ## simulation
  for k = 1 : trows
    y(k, :) = C * x  +  D * u(k, :)';
    x_arr(k, :) = x;
    x = F * x  +  G * u(k, :)';
  endfor

  if (! nargout)  # plot information
    str = "Linear Simulation Results";
    outname = get (sys, "outname");

    if (isempty (outname) || isequal ("", outname{:}))
      outname = strseq ("y_", 1 : length (outname));
    else
      outname = __markemptynames__ (outname);
    endif

    if (digital)  # discrete system
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
    else  # continuous system
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
  else  # return values
    y_r = y;
    t_r = t;
    x_r = x_arr;
  endif

endfunction
