## Copyright (C) 1996, 1998, 2000, 2003, 2004, 2005, 2006, 2007
##               Auburn University.  All rights reserved.
## Copyright (C) 2009 Lukas Reichlin.  All rights reserved.
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn{Function File} {[@var{y}, @var{t}, @var{x}] =} initial (@var{sys}, @var{x0})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} initial (@var{sys}, @var{x0}, @var{tfinal})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} initial (@var{sys}, @var{x0}, @var{tfinal}, @var{dt})
## Initial condition response of state-space model.
## If no output arguments are given, the response is printed on the screen;
## otherwise, the response is computed and returned.
##
## @strong{Inputs}
## @table @var
## @item sys
## System data structure. Must be either purely continuous or discrete;
## see @code{is_digital}.
## @item x0
## Vector of initial conditions for each state.
## @item tfinal
## Optional simulation horizon. If not specified, it will be calculated by
## the poles of the system to reflect adequately the response transients.
## @item dt
## Optional sampling time. Be sure to choose it small enough to capture transient
## phenomena. If not specified, it will be calculated by the poles of the system.
## @end table
##
## @strong{Outputs}
## @table @var
## @item y
## Output response array. Has as many rows as time samples (length of t)
## and as many columns as outputs.
## @item t
## Time row vector.
## @item x
## State trajectories array. Has length(t) rows and as many columns as states.
## @end table
##
## @seealso{impulse, lsim, step}
## @example
## @group
##                    .
## Continuous Time:   x = A x ,   y = C x ,   x(0) = x0
##
## Discrete Time:   x[k+1] = A x[k] ,   y[k] = C x[k] ,   x[0] = x0
## @end group
## @end example
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@swissonline.ch>
## Created: August 16, 2009
## based on __stepimp__.m of Kai P. Mueller and A. Scottedward Hodel
## Version: 0.1

function [y_r, t_r, x_r] = initial (sys, x_0, t_final, dt)

  ## check whether arguments are OK
  if (nargin < 2 || nargin > 4)
    print_usage ();
  endif

  if (! isstruct (sys))
    error ("initial: first argument must be a system data structure");
  endif

  if (! isvector (x_0))
    error ("initial: second argument must be a vector");
  endif

  ## get system information
  sys = sysupdate (sys, "ss");
  digital = is_digital (sys, 2);
  [n_c, n_d, n_in, n_out] = sysdimensions (sys);
  n_st = n_c + n_d;  # number of states

  if (digital == -1)
    error ("initial: system must be either purely continuous or purely discrete");
  endif

  if (n_st != length (x_0))
    error ("initial: x0 must be a vector with %d elements", n_st);
  endif

  ## code adapted from __stepimp__.m
  if (nargin < 3)
    ## we have to compute the time when the system reaches steady state
    ## and the step size
    eigw = eig (sys2ss (sys));

    if (digital)
      t_sam = sysgettsam (sys);

      ## perform bilinear transformation on poles in z
      for k = 1 : n_d
        pole = eigw(k);
        if (abs (pole + 1) < 1.0e-10)
          eigw(k) = 0;
        else
          eigw(k) = 2 / t_sam * (pole - 1) / (pole + 1);
        endif
      endfor
    endif

    ## remove poles near zero from eigenvalue array eigw
    nk = n_st;
    for k = 1 : n_st
      if (abs (real (eigw(k))) < 1.0e-10)
        eigw(k) = 0;
        nk = nk - 1;
      endif
    endfor

    if (nk == 0)
      t_final = 10;
      dt = t_final / 1000;
    else
      eigw = eigw(find (eigw));
      eigw_max = max (abs (eigw));
      dt = 0.2 * pi / eigw_max;

      eigw_min = min (abs (real (eigw)));
      t_final = 5.0 / eigw_min;

      ## round up
      yy = 10^(ceil (log10 (t_final)) - 1);
      t_final = yy * ceil (t_final / yy);

      n = t_final / dt;

      if (n < 50)
        dt = t_final / 50;
      endif

      if (n > 2000)
        dt = t_final / 2000;
      endif
    endif
  endif
  ## end of adapted code

  if (digital)
    dt = sysgettsam (sys);

    if (nargin == 4)
      warning ("initial: fourth argument has no effect on sampling time of digital system");
    endif
  else
    sys = c2d (sys, dt);
  endif

  t = (0 : dt : t_final)';
  l_t = length (t);

  [F, G, C, D] = sys2ss (sys);

  ## preallocate memory
  y = zeros (l_t, n_out);
  x_vec = zeros (l_t, n_st);

  ## make sure that x is a row vector
  x = reshape (x_0, length (x_0), 1);

  ## simulation
  for k = 1 : l_t
    y(k, :) = C * x;
    x_vec(k, :) = x;
    x = F * x;
  endfor

  if (nargout == 0)  # plot information
    if (digital)  # discrete system
      for k = 1 : n_out
        subplot (n_out, 1, k)
        stairs (t, y(:, k))
        grid on
        if (k == 1)
          title ("Response to Initial Conditions")
        endif
        ylabel (sprintf ("Amplitude %s", sysgetsignals (sys, "out", k, 1)))
      endfor
      xlabel ("Time [s]")
    else  # continuous system
      for k = 1 : n_out
        subplot (n_out, 1, k)
        plot (t, y(:, k))
        grid on
        if (k == 1)
          title ("Response to Initial Conditions")
        endif
        ylabel (sprintf ("Amplitude %s", sysgetsignals (sys, "out", k, 1)))
      endfor
      xlabel ("Time [s]")
    endif  
  else  # return values
    y_r = y;
    t_r = t;
    x_r = x_vec;
  endif

endfunction


## TODO: Add a purely continuous system test
## TODO: Add a purely discrete system test
