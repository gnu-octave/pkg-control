## Copyright (C) 2024 Fabio Di Iorio <diiorio.fabio@gmail.com>
##
## This file is part of the Control package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not,
## see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {@var{reg} =} lqg (@var{sys}, @var{QXU}, @var{QWV})
## @deftypefnx {Function File} {@var{reg} =} lqg (@var{sys}, @var{QXU}, @var{QWV}, @var{QI})
## @deftypefnx {Function File} {@var{reg} =} lqg (@var{sys}, @var{QXU}, @var{QWV}, @var{QI}, @var{sensors}, @var{known})
## @deftypefnx {Function File} {@var{reg} =} lqg (@var{sys}, @var{QXU}, @var{QWV}, @var{}, @var{sensors}, @var{known})
## Linear-quadratic gaussian (LQG) design
##
## @strong{Inputs}
## @table @var
## @item sys
## Continuous or discrete-time @acronym{LTI} model (m inputs, n states, p outputs).
## @item QXU
## State and input weighting matrix (n+m-by-n+m).
## @item QWV
## Process and measurement noise covariance matrix (n-by-n).
## @item QI
## Optional output weighting matrix for LQG servo control with integral action (p-by-p).  If @var{QI} is not specified, the LQG regulator is computed
## @end table
##
## @strong{Outputs}
## @table @var
## @item reg
## LQG regulator or controller as dynamic compensator. Connect with positive feedback.
## @end table
##
## @strong{Equations}
## @example
## @group
## .
## x = A x + B u,   x(0) = x0
##
##                   1   T
## J(x0) = E{ lim   --- INT ([x', u'] Qxu [x u]' + xi' Qi xi) dt}
##           T->inf  T   0
##
## @end group
## @end example
## @seealso{lqr, kalman, lqi}
## @end deftypefn

## Author: Fabio Di Iorio <diiorio.fabio@gmail.com>
## Created: July 2024
## Version: 0.1

function [reg] = lqg (sys, QXU, QWV, QI = [])

  if (nargin < 3 || nargin > 4)
    print_usage ();
  endif

  if (isa (sys, "lti"))
    [a, b, c, d, e, Ts] = dssdata (sys, []);
  else
    print_usage ();
  endif

  if (~issquare (QXU) || ~issquare (QWV) || (~isempty(QI) && ~issquare(QI)))
    print_usage ();
  endif

  [n, m] = size(b);
  [p, ~] = size(c);

  % LQR data
  Q = QXU(1:n, 1:n);
  R = QXU(n+1:end, n+1:end);
  S = QXU(1:n, n+1:end);

  % observer data
  W = QWV(1:n, 1:n);
  V = QWV(n+1:end, n+1:end);
  N = QWV(1:n, n+1:end);

  if (isempty(QI))
    K = lqr(sys, Q, R, S);
  else
    K = lqi(sys, blkdiag(Q, QI), R, vertcat(S,zeros(p,m)));
  endif

  % add noise to system as additional inputs, to use with "kalman" function
  Bn = [b eye(n)];
  sys_noisy = ss(a, Bn, c, d, Ts);

  [est, L, ~] = kalman(sys_noisy, W, V, N, 1:p, 1:m);

  % regulator case
  if isempty(QI)
    reg = ss(a-b*K-L*c-L*d*K, L, -K, 0,Ts);

    % set variables names
    [inn, stn, outn, ing, outg] = get (sys, "inname", "stname", "outname", "ingroup", "outgroup");
    stname = __labels__ (stn, "xhat");
    outname = vertcat (__labels__ (outn(1:m), "u"));
    inname = vertcat (__labels__ (outn(1:p), "y"));
    reg = set (reg, "inname", inname, "stname", stname, "outname", outname);
  else
    % servo case
    if isct(sys)
      reg = ss([a-b*K(:,1:(end-p))-L*c+L*d*K(:,1:(end-p)) -b*K(:,(n+1):end)+L*d*K(:,(n+1):end); zeros(p,n) zeros(p,p)], [zeros(n,p) L; ones(p,1) -1.*ones(p,1)], -K, 0,Ts);
    else
      reg = ss([a-b*K(:,1:(end-p))-L*c+L*d*K(:,1:(end-p)) -b*K(:,(n+1):end)+L*d*K(:,(n+1):end); zeros(p,n) eye(p,p)], [zeros(n,p) L; Ts.*ones(p,1) -Ts.*ones(p,1)], -K, 0,Ts);
    endif

    % set variables names
    [inn, stn, outn, ing, outg] = get (sys, "inname", "stname", "outname", "ingroup", "outgroup");
    stname = vertcat (__labels__ (stn(1:n), "xhat"), __labels__ (outn(1:p), "xi"));
    outname = vertcat (__labels__ (outn(1:m), "u"));
    inname = vertcat (__labels__ (outn(1:p), "r"),__labels__ (outn(1:p), "y"));
    reg = set (reg, "inname", inname, "stname", stname, "outname", outname);
  endif

endfunction

%!test
%! G = zpk([], [-10 -1 -100], 2000);
%! sys = ss(G);
%! [A B C D] = ssdata(sys);
%! Q = eye(3);
%! QI = 100;
%! QXU = blkdiag (Q, 1);
%! QWV = eye(4);
%! reg = lqg(sys, QXU, QWV);
%! assert(real(eig(feedback(reg, sys, 1)))<0);
%! reg=lqg(sys,QXU,QWV,QI);
%! assert(real(eig(feedback(reg, sys, 2, 1, 1)))<0);
%! Ts = 0.01;
%! Gz = zpk([], [-0.1 0.05 0.004], 3, Ts);
%! sysz = ss(Gz);
%! regz = lqg(sysz, QXU, QWV);
%! assert(abs(eig(feedback(regz, sysz, 1)))<1);
%! regz=lqg(sysz, QXU, QWV, QI);
%! assert(abs(eig(feedback(regz, sysz, 2, 1, 1)))<1);
