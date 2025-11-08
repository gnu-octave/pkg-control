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
## @deftypefn {Function File} {@var{reg} =} reg (@var{sys}, @var{k}, @var{l})
## Form regulator from state-feedback and estimator gains
##
## @strong{Inputs}
## @table @var
## @item sys
## State-space model of the plant
## @item k
## State-feedback gain
## @item l
## Estimator gain
## @end table
##
## @strong{Outputs}
## @table @var
## @item reg
## Dynamic compensator. Connect with positive feedback.
## @end table
##
## @strong{Equations}
## @seealso{place}
## @end deftypefn

## Author: Fabio Di Iorio <diiorio.fabio@gmail.com>
## Created: July 2024
## Version: 0.1

function regsys = reg (sys, k, l)

  ## TODO: implement variant with additional known inputs:

    [m, n] = size(k);
    [~, p] = size(l);

    if (nargin  ~= 3)
      print_usage ();
    endif

    if (~isnumeric(k) || ~isnumeric(l))
      print_usage ();
    endif

    if (n<1 || m<1 || p<1)
      print_usage ();
    endif

    if (~isa (sys, "lti"))
      print_usage ();
    endif

    A = sys.a;
    C = sys.c;
    B = sys.b;
    D = sys.d;
    Ts = sys.ts;
    regsys = ss(A-B*k-l*C-l*D*k, l, -k, 0,Ts);
    % set variables names
    [inn, stn, outn, ing, outg] = get (sys, "inname", "stname", "outname", "ingroup", "outgroup");
    stname = __labels__ (stn, "xhat");
    outname = cell(m,1);
    for i=1:m
      outname{i,1} = strcat("u",num2str(i));
    endfor
    inname = cell(p,1);
    for i=1:p
      inname{i,1} = strcat("y",num2str(i));
    endfor
    regsys = set (regsys, "inname", inname, "stname", stname, "outname", outname);

endfunction

%!test
%! G = zpk([],[-10 -1 -100], 2000);
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
%! Creg = reg(sys, K, L1);
%! assert(real(eig(feedback(Creg, sys, 1)))<0);
%! Ts = 0.01;
%! Gz = zpk([],[-0.1 0.05 0.004], 3, Ts);
%! sysz = ss(Gz);
%! kz = lqr(sysz, Q, R, S);
%! Bn = [sysz.b eye(n)];
%! sys_noisyz = ss(sysz.a, Bn, sysz.c, sysz.d, sysz.ts);
%! [estz, L1z, ~] = kalman(sys_noisyz, W, V, N, 1:p, 1:m);
%! Cz = reg(sysz, kz, L1z);
%! assert(abs(eig(feedback(Cz,sysz,1)))<1);
