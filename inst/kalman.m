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
## @deftypefn {Function File} {[@var{est}, @var{g}, @var{x}] =} kalman (@var{sys}, @var{Q}, @var{R})
## @deftypefnx {Function File} {[@var{est}, @var{g}, @var{x}] =} kalman (@var{sys}, @var{Q}, @var{R}, @var{S})
## @deftypefnx {Function File} {[@var{est}, @var{g}, @var{x}] =} kalman (@var{sys}, @var{Q}, @var{R}, @var{[]}, @var{sensors}, @var{known})
## @deftypefnx {Function File} {[@var{est}, @var{g}, @var{x}] =} kalman (@var{sys}, @var{Q}, @var{R}, @var{S}, @var{sensors}, @var{known})
## @deftypefnx {Function File} {[@var{est}, @var{g}, @var{x}] =} kalman (@var{sys}, @var{Q}, @var{R}, @var{[]}, @var{sensors}, @var{known}, @var{type})
## @deftypefnx {Function File} {[@var{est}, @var{g}, @var{x}] =} kalman (@var{sys}, @var{Q}, @var{R}, @var{S}, @var{sensors}, @var{known}, @var{type})
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
## Optional cross term covariance. Default value is 0.
## @item sensors
## Indices of measured output signals y from @var{sys}. If omitted or empty, all outputs are measured.
## @item known
## Indices of known input signals u (deterministic) to @var{sys}.  All other inputs to @var{sys}
## are assumed stochastic. If argument @var{known} is omitted or empty, the first m-l inputs to @var{sys}
## are known, where m is the total number of inputs to @var{sys} and l is the size of the quadratic
## matrix @var{Q}.
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

function [est, K, X] = kalman (sys, Q, R, varargin)

  if (nargin < 3 || nargin > 7 || ! isa (sys, "lti"))
    print_usage ();
  endif

  ## System in state space
  [A, B, C, D, E] = dssdata (sys, []);

  ## Optional parameters:
  ## The new default for known inputs (deterministic) are not
  ## backward compatible because previously, all inputs are
  ## assumed to be stochastic if is empty. However, this would
  ## result in an error if the number if the number of stochastic
  ## inputs would not match the dimension of Q. For all valid cases
  ## the new version is compatibel to the old one.
  S = [];
  sensors = 1 : rows (C);
  deterministic = 1 : columns (B) - size (Q,1);   # first inputs deterministic
  type = 'delayed';

  varidxoff = 3;        # offset no. of variable and fixed input arguments
  argidx = varidxoff;   # index of last input argument

  if (nargin > argidx++)
    S = varargin{argidx-varidxoff};
    if (nargin > argidx++)
      if (! isempty (varargin{argidx-varidxoff}))
        sensors = varargin{argidx-varidxoff};
      endif
      if (nargin > argidx++)
        if (! isempty (varargin{argidx-varidxoff}))
          deterministic = varargin{argidx-varidxoff};
        endif
        if (nargin > argidx++)
          if (isct (sys))
            warning ("kalman: ignoring 'type' parameter for continuous-time estimator\n");
          else
            type = varargin{argidx-varidxoff};
          endif
        endif
      endif
    endif
  endif

  m   = columns (B);
  m_d = length (deterministic);
  m_s = size (Q,1);
  p = length (sensors);
  p_s = size (R,1);

  ## plausibility check for parameters
  if (! issquare (Q))
    error ("kalman: second argument Q must be square\n");
  endif
  if (! issquare (R))
    error ("kalman: third argument R must be square\n");
  endif
  if (m_s != m - m_d)
    error ("kalman: number of stochastic inputs (%d) does not match size %d of Q\n",...
           m - m_d, m_s);
  endif
  if (p_s != p)
    error ("kalman: size %d of measurment noise does not match size %d of R\n",...
           p, p_s);
  endif
  if ((! isempty (S)) && (size (S) != [m_s,p_s]))
    error ("kalman: size [%d,%d] of S does not match size %d of Q and %d of R\n",...
           size(S,1), size(S,2), m_s, p_s);
  endif

  ## matrices for Kalman filter design
  stochastic = setdiff (1 : columns (B), deterministic);

  C = C(sensors, :);
  G = B(:, stochastic);
  H = D(sensors, stochastic);

  if (isempty (S))
    Rbar = R + H*Q*H.';
    Sbar = G * Q*H.';
  else
    Rbar = R + H*Q*H.'+ H*S + S.'*H.';
    Sbar = G * (Q*H.' + S);
  endif

  if (isct (sys))
    [X, L, K] = care (A.', C.', G*Q*G.', Rbar, Sbar, E.');
  else
    [X, L, K] = dare (A.', C.', G*Q*G.', Rbar, Sbar, E.');
  endif

  K = K.';
  est = estim (sys, K, sensors, deterministic, type);

endfunction

%!test
%! sys = ss (-2, 1, 1, 3);
%! [est, g, x] = kalman (sys, 1, 1, 1);
%! [a, b, c, d] = ssdata (est);
%! m = [a, b; c, d];
%! m_exp = [-2.25, 0.25; 1, 0; 1, 0];
%! g_exp = 0.25;
%! x_exp = 0;
%! assert (m, m_exp, 1e-2);
%! assert (g, g_exp, 1e-2);
%! assert (x, x_exp, 1e-2);

%!shared n, nw, A, B, C, D, Bw, Dw, sys, Q, R, N, Pinf, K
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
%! D = [1];
%! Dw = zeros(1,nw);
%! D = [D Dw];
%! sys = ss(A,B,C,D,1);
%! Q = eye(nw,nw);
%! R = 1;
%! N = [];
%! PP = eye(3,3);
%! # asymptotic filter equations
%! for i = 1:10
%!   Pinf = A*PP*A' + Bw*Q*Bw';
%!   K  = Pinf*C'*inv(C*Pinf*C'+R);
%!   PP =  (eye(3,3)-K*C)*Pinf;
%! endfor;

%!test
%! [kalmf,L,P] = kalman(sys,Q,R,N,1,1);
%! assert (Pinf, P, 4e-4);
%! assert (inv(A)*L, K, 3e-4);
%! assert (kalmf.a, A-L*C, 1e-4);

%!test
%! [kalmfc,Lc,Pc] = kalman(sys,Q,R,N,1,1,'current');
%! assert (Pinf, Pc, 4e-4);
%! assert (inv(A)*Lc, K, 3e-4);
%! assert (kalmfc.a, A-A*K*C, 1e-4);

