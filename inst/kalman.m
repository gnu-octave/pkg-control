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
## @deftypefn {Function File} {[@var{est}, @var{g}, @var{x}] =} kalman (@var{sys}, @var{q}, @var{r})
## @deftypefnx {Function File} {[@var{est}, @var{g}, @var{x}] =} kalman (@var{sys}, @var{q}, @var{r}, @var{s})
## @deftypefnx {Function File} {[@var{est}, @var{g}, @var{x}] =} kalman (@var{sys}, @var{q}, @var{r}, @var{[]}, @var{sensors}, @var{known})
## @deftypefnx {Function File} {[@var{est}, @var{g}, @var{x}] =} kalman (@var{sys}, @var{q}, @var{r}, @var{s}, @var{sensors}, @var{known})
## @deftypefnx {Function File} {[@var{est}, @var{g}, @var{x}] =} kalman (@var{sys}, @var{q}, @var{r}, @var{[]}, @var{sensors}, @var{known}, @var{type})
## @deftypefnx {Function File} {[@var{est}, @var{g}, @var{x}] =} kalman (@var{sys}, @var{q}, @var{r}, @var{s}, @var{sensors}, @var{known}, @var{type})
## Design Kalman estimator for @acronym{LTI} systems.
##
## @strong{Inputs}
## @table @var
## @item sys
## Nominal plant model.
## @item q
## Covariance of white process noise.
## @item r
## Covariance of white measurement noise.
## @item s
## Optional cross term covariance.  Default value is 0.
## @item sensors
## Indices of measured output signals y from @var{sys}.  If omitted, all outputs are measured.
## @item known
## Indices of known input signals u (deterministic) to @var{sys}.  All other inputs to @var{sys}
## are assumed stochastic.  If argument @var{known} is omitted, no inputs u are known.
## @item type
## Type of the estimator for discrete-time systems. If set to 'delayed' the current
## estimation is based on y(k-1), if set to 'current' the current estimation is
## based on the lates mesaruement y(k). If omitted, the 'delayed' version is created.
## @end table
##
## @strong{Outputs}
## @table @var
## @item est
## State-space model of the Kalman estimator.
## @item g
## Estimator gain.
## @item x
## Solution of the Riccati equation.
## @end table
##
## @strong{Block Diagram}
## @example
## @group
##                                  u  +-------+         ^
##       +---------------------------->|       |-------> y
##       |    +-------+     +       y  |  est  |         ^
## u ----+--->|       |----->(+)------>|       |-------> x
##            |  sys  |       ^ +      +-------+
## w -------->|       |       |
##            +-------+       | v
##
## Q = cov (w, w')     R = cov (v, v')     S = cov (w, v')
## @end group
## @end example
##
## @seealso{care, dare, estim, lqr}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: November 2009
## Version: 0.3

function [est, k, x] = kalman (sys, q, r, s = [], sensors = [], deterministic = [], varargin)

  if (nargin < 3 || nargin > 7 || ! isa (sys, "lti"))
    print_usage ();
  endif

  ## optional parameters
  type = 'delayed';
  varidx = 0;

  if (nargin > 6)
    if (isct (sys))
      warning ("kalman: ignoring 'type' parameter for continuous-time estimator\n");
    else
      type = varargin{++varidx};
    endif
  endif

  [a, b, c, d, e] = dssdata (sys, []);

  if (isempty (sensors))
    sensors = 1 : rows (c);
  endif

  stochastic = setdiff (1 : columns (b), deterministic);

  c = c(sensors, :);
  g = b(:, stochastic);
  h = d(sensors, stochastic);

  if (isempty (s))
    rbar = r + h*q*h.';
    sbar = g * q*h.';
  else
    rbar = r + h*q*h.'+ h*s + s.'*h.';
    sbar = g * (q*h.' + s);
  endif

  if (isct (sys))
    [x, l, k] = care (a.', c.', g*q*g.', rbar, sbar, e.');
  else
    [x, l, k] = dare (a.', c.', g*q*g.', rbar, sbar, e.');
  endif

  k = k.';
  est = estim (sys, k, sensors, deterministic, type);

endfunction

%!test
%!shared m, m_exp, g, g_exp, x, x_exp
%! sys = ss (-2, 1, 1, 3);
%! [est, g, x] = kalman (sys, 1, 1, 1);
%! [a, b, c, d] = ssdata (est);
%! m = [a, b; c, d];
%! m_exp = [-2.25, 0.25; 1, 0; 1, 0];
%! g_exp = 0.25;
%! x_exp = 0;
%!assert (m, m_exp, 1e-2);
%!assert (g, g_exp, 1e-2);
%!assert (x, x_exp, 1e-2);

%!test
%!shared P, Pc, Pinf, L, Lc, K, A, C, kalmf, kalmfc
%! n = 3;
%! nw = 2;
%!
%! A = [1.1269   -0.4940    0.1129
%!      1.0000         0         0
%!           0    1.0000         0];
%! B = [-0.3832
%!       0.5919
%!       0.5191];
%! Bw = rand(n,nw);
%! B = [B Bw];
%! C = [1 0 0];
%! D = [0];
%! Dw = zeros(1,nw);
%! sys = ss(A,B,C,D,1);
%! Q = eye(nw,nw);
%! R = 1;
%! N = [];
%! [kalmf,L,P] = kalman(sys,Q,R,N,1,1);
%! [kalmfc,Lc,Pc] = kalman(sys,Q,R,N,1,1,'current');
%! PP = eye(3,3);
%! # asymptotic filter equations
%! for i = 1:10
%!   Pinf = A*PP*A' + Bw*Q*Bw';
%!   K  = Pinf*C'*inv(C*Pinf*C'+R);
%!   PP =  (eye(3,3)-K*C)*Pinf;
%! end;
%!assert (Pinf, P, 1e-4);
%!assert (Pinf, Pc, 1e-4);
%!assert (inv(A)*L, K, 1e-4);
%!assert (inv(A)*Lc, K, 1e-4);
%!assert (kalmf.a, A-L*C, 1e-4);
%!assert (kalmfc.a, A-A*K*C, 1e-4);

