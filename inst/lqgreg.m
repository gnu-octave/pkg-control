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
## @deftypefn {Function File} {@var{reg} =} lqgreg (@var{kest}, @var{k})
## Form LQG regulator
##
## @strong{Inputs}
## @table @var
## @item kest
## Kalman estimator
## @item k
## State-feedback gain
## @end table
##
## @strong{Outputs}
## @table @var
## @item reg
## LQG regulator as dynamic compensator. Connect with positive feedback.
## @end table
##
## @strong{Equations}
## @seealso{lqr, kalman, lqg}
## @end deftypefn

## Author: Fabio Di Iorio <diiorio.fabio@gmail.com>
## Created: July 2024
## Version: 0.1

function [reg] = lqgreg (kest, k)

  ## TODO: implement variant with additional known inputs:
  ## reg = lqgreg (kest, k, controls)

    if (isa (kest, "lti"))
      [a, b, c, d, e, Ts] = dssdata (kest, []);
    else
      print_usage ();
    endif

    [m, ~] = size(k);

    L = kest.b(:, m+1:end);
    [~, p] = size(L);
    C = kest.c(1:p, :);
    B = kest.b(:, 1:m);
    D = kest.d(1:p, 1:m);
    reg = ss(a-B*k-L*C-L*D*k, L, -k, 0,Ts);
    % set variables names
    [inn, stn, outn, ing, outg] = get (kest, "inname", "stname", "outname", "ingroup", "outgroup");
    stname = __labels__ (stn, "xhat");
    outname = cell(m,1);
    for i=1:m
      outname{i,1} = strcat("u",num2str(i));
    endfor
    inname = cell(p,1);
    for i=1:p
      inname{i,1} = strcat("y",num2str(i));
    endfor
    reg = set (reg, "inname", inname, "stname", stname, "outname", outname);

endfunction

%!test
%! G=zpk([], [-10 -1 -100], 2000);
%! sys = ss(G);
%! [n, m] = size(sys.b);
%! [p, ~] = size(sys.c);
%! Q = eye(3);
%! R = 1;
%! S = zeros(3, 1);
%! W = eye(3);
%! V = 1;
%! N = zeros(3, 1);
%! K = lqr(sys, Q, R, S);
%! Bn = [sys.b eye(n)];
%! sys_noisy = ss(sys.a, Bn, sys.c, sys.d, sys.ts);
%! [est, L1, ~] = kalman(sys_noisy, W, V, N, 1:p, 1:m);
%! reg = lqgreg(est,K);
%! assert(real(eig(feedback(reg, sys, 1))) < 0);
%! Ts = 0.01;
%! Gz=zpk([], [-0.1 0.05 0.004], 3, Ts);
%! sysz = ss(Gz);
%! kz = lqr(sysz, Q, R, S);
%! Bn = [sysz.b eye(n)];
%! sys_noisyz = ss(sysz.a, Bn, sysz.c, sysz.d, sysz.ts);
%! [estz, L1, ~] = kalman(sys_noisyz, W, V, N, 1:p, 1:m);
%! regz = lqgreg(estz, kz);
%! assert(abs(eig(feedback(regz, sysz, 1))) < 1);
