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
## @deftypefn {Function File} {@var{reg} =} lqgtrack (@var{kest}, @var{k})
## Form LQG servo controller
##
## @strong{Inputs}
## @table @var
## @item kest
## Kalman estimator
## @item k
## State-feedback gain, including integrator states (m x n+p)
## @end table
##
## @strong{Outputs}
## @table @var
## @item reg
## LQG servo controller. Connect with positive feedback.
## @end table
##
## @strong{Equations}
## @seealso{lqr, kalman, lqg, lqgreg}
## @end deftypefn

## Author: Fabio Di Iorio <diiorio.fabio@gmail.com>
## Created: July 2024
## Version: 0.1

function [reg] = lqgtrack (kest, k)

  %% TODO: implement variant with additional known inputs:
  %% reg = lqgtrack (kest, k, controls)
  %  TODO: implement variant with 1dof
  %% reg = lqgtrack (kest, k, '1dof')

    if (isa (kest, "lti"))
      [a, b, c, d, e, Ts] = dssdata (kest, []);
    else
      print_usage ();
    endif

    [m, ~] = size(k);

    L = kest.b(:, m+1:end);
    [n, p] = size(L);
    C = kest.c(1:p, :);
    B = kest.b(:, 1:m);
    D = kest.d(1:p, 1:m);
    if isct(kest)
      reg = ss([a-B*k(:,1:(end-p))-L*C+L*D*k(:,1:(end-p)) -B*k(:,(n+1):end)+L*D*k(:,(n+1):end); zeros(p,n) zeros(p,p)], [zeros(n,p) L; ones(p,1) -1.*ones(p,1)], -k, 0,Ts);
    else
      reg = ss([a-B*k(:,1:(end-p))-L*C+L*D*k(:,1:(end-p)) -B*k(:,(n+1):end)+L*D*k(:,(n+1):end); zeros(p,n) eye(p,p)], [zeros(n,p) L; Ts.*ones(p,1) -Ts.*ones(p,1)], -k, 0,Ts);
    endif
    % set variables names
    [inn, stn, outn, ing, outg] = get (kest, "inname", "stname", "outname", "ingroup", "outgroup");
    stname = cell(n+p,1);
    for i=1:n
      stname{i,1} = strcat("xhat",num2str(i));
    endfor
    for i=1:p
      stname{n+i,1} = strcat("xi",num2str(i));
    endfor
    outname = cell(m,1);
    for i=1:m
      outname{i,1} = strcat("u",num2str(i));
    endfor
    inname = cell(2*p,1);
    for i=1:p
      inname{i,1} = strcat("r",num2str(i));
    endfor
    for i=1:p
      inname{p+i,1} = strcat("y",num2str(i));
    endfor
    reg = set (reg, "inname", inname, "stname", stname, "outname", outname);

endfunction

%!test
%! G=zpk([],[-10 -1 -100],2000);
%! sys=ss(G);
%! [n, m] = size(sys.b);
%! [p, ~] = size(sys.c);
%! Q=eye(3+p);
%! R = 1;
%! S = zeros(3+p,1);
%! W = eye(3);
%! V = 1;
%! N = zeros(3,1);
%! K = lqi(sys, Q, R, S);
%! Bn = [sys.b eye(n)];
%! sys_noisy = ss(sys.a, Bn, sys.c, sys.d, sys.ts);
%! [est, L1, ~] = kalman(sys_noisy, W, V, N, 1:p, 1:m);
%! reg=lqgtrack(est,K);
%! assert(real(eig(feedback(reg,sys,2,1,1)))<0);
%! Ts=0.01;
%! Gz=zpk([],[-0.1 0.05 0.004],3,Ts);
%! sysz=ss(Gz);
%! kz = lqi(sysz, Q, R, S);
%! Bn = [sysz.b eye(n)];
%! sys_noisyz = ss(sysz.a, Bn, sysz.c, sysz.d, sysz.ts);
%! [estz, L1, ~] = kalman(sys_noisyz, W, V, N, 1:p, 1:m);
%! regz=lqgtrack(estz,kz);
%! assert(abs(eig(feedback(regz,sysz,2,1,1)))<1);
