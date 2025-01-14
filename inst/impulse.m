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
## @deftypefn{Function File} {} impulse (@var{sys})
## @deftypefnx{Function File} {} impulse (@var{sys1}, @var{sys2}, @dots{}, @var{sysN})
## @deftypefnx{Function File} {} impulse (@var{sys1}, @var{'style1'}, @dots{}, @var{sysN}, @var{'styleN'})
## @deftypefnx{Function File} {} impulse (@var{sys1}, @dots{}, @var{t})
## @deftypefnx{Function File} {} impulse (@var{sys1}, @dots{}, @var{tfinal})
## @deftypefnx{Function File} {} impulse (@var{sys1}, @dots{}, @var{tfinal}, @var{dt})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} impulse (@var{sys})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} impulse (@var{sys}, @var{t})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} impulse (@var{sys}, @var{tfinal})
## @deftypefnx{Function File} {[@var{y}, @var{t}, @var{x}] =} impulse (@var{sys}, @var{tfinal}, @var{dt})
## Impulse response of @acronym{LTI} system.
## If no output arguments are given, the response is printed on the screen.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} model.
## @item t
## Time vector.  Should be evenly spaced.  If not specified, it is calculated by
## the poles of the system to reflect adequately the response transients.
## @item tfinal
## Optional simulation horizon.  If not specified, it is calculated by
## the poles of the system to reflect adequately the response transients.
## @item dt
## Optional sampling time.  Be sure to choose it small enough to capture transient
## phenomena.  If not specified, it is calculated by the poles of the system.
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
## Time row vector.
## @item x
## State trajectories array.  Has @code{length (t)} rows and as many columns as states.
## @end table
##
## @strong{Remark}
##
## For the impulse response of a discrete-time system, the input
## sequence @{1/T,0,0,0,...@} and not the unit impulse is considered.
##
## @seealso{initial, lsim, step}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 1.0

function [y_r, t_r, x_r] = impulse (varargin)

  if (nargin == 0)
    print_usage ();
  endif

  names = cell (1,nargin);
  for i = 1:nargin
    names{i} = inputname (i);
  end

  [y, t, x] = __time_response__ ("impulse", varargin, names, nargout);

  if (nargout)
    y_r = y{1};
    t_r = t{1};
    x_r = x{1};
  endif

endfunction


## Test the analogue response to an impulse input.
## In this case the system in converted into state space
## and then discretized. In the test below, the system is
## discretized as transfer function and then transferred into
## state space for simulating the impulse response. The results
## differ in a quantity larger than 2*eps, therefore 10*eps is
## chosen as tolerance.
%!test
%! t=0:1:4;
%! sys=tf(1,[1 2 2 1]);
%! y=impulse(sys,t);
%! assert (y(1), 0, eps);
%! assert (y(2), 0.241686482894434, 2*eps);
%! assert (y(3), 0.404040547757057, 2*eps);
%! assert (y(4), 0.307384479794317, 2*eps);
%! assert (y(5), 0.121908527560869, 6*eps);  # error propagation

## Test the discrete response to an impulse input.
%!test
%! t=0:1:4;
%! sys=tf(1,[1 2 2 1]);
%! sys2=c2d(sys, 1, "impulse");
%! y=impulse(sys2,t);
%! assert (y(1), 0, eps);
%! assert (y(2), 0.241686482894434, 2*eps);
%! assert (y(3), 0.404040547757057, 3*eps);
%! assert (y(4), 0.307384479794317, 7*eps);
%! assert (y(5), 0.121908527560869, 7*eps);  # error propagation

## test from bug
%!test
%! s = tf("s");
%! R = 1/s;
%! y= impulse(R);
%! assert (y(1), 1, eps);

%!demo
%! clf;
%! s = tf('s');
%! g = 1/(2*s^2+3*s+4);
%! impulse(g);
%! title ("Impulse response of a PT2 transfer function");

%!demo
%! clf;
%! s = tf('s');
%! g = 1/(2*s^2+3*s+4);
%! h = c2d(g,0.1);
%! impulse(h);
%! title ("Impulse response of a discretized PT2 transfer function");
