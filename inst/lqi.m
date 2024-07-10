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
## @deftypefn {Function File} {[@var{g}, @var{x}, @var{l}] =} lqi (@var{sys}, @var{q}, @var{r})
## @deftypefnx {Function File} {[@var{g}, @var{x}, @var{l}] =} lqr (@var{sys}, @var{q}, @var{r}, @var{s})
## Linear-quadratic integral control.
##
## @strong{Inputs}
## @table @var
## @item sys
## Continuous or discrete-time @acronym{LTI} model (m inputs, n states, p outputs).
## @item q
## State weighting matrix (n+p-by-n+p).
## @item r
## Input weighting matrix (m-by-m).
## @item s
## Optional cross term matrix (n+p-by-m).  If @var{s} is not specified, a zero matrix is assumed.
## @end table
##
## @strong{Outputs}
## @table @var
## @item g
## State feedback matrix (m-by-n).
## @item x
## Unique stabilizing solution of the continuous-time Riccati equation (n+p-by-n+p).
## @item l
## Closed-loop poles (n-by-1).
## @end table
##
## @strong{Equations}
## @example
## @group
## .
## x = A x + B u,   x(0) = x0
##
##         inf
## J(x0) = INT (z' Q z  +  u' R u  +  2 z' S u)  dt
##          0
##
## z = [x; xi]
#  z is the augmented state and xi the integral of the error e = r - y
## L = eig (A - B*G)
## @end group
## @end example
## @seealso{care, dare, dlqr}
## @end deftypefn

## Author: Fabio Di Iorio <diiorio.fabio@gmail.com>
## Created: June 2024
## Version: 0.1

function [g, x, l] = lqi (sys, q, r, s = [])

  if (nargin < 3 || nargin > 4)
    print_usage ();
  endif

  if ~isa(sys, 'ss')
    print_usage ();
  endif

[A B C D Ts] = ssdata(sys);

[n, m] = size(B);
[p, ~] = size(C);

% create augmented plant with integrator states

if Ts==0
  % continuous time case
  Ag = [A zeros(n,p); -C zeros(p,p)];
  Bg = [B; -D];
  Cg = [C zeros(p, p)];
  sysg = ss(Ag, Bg, Cg, D);
else
  % discrete time case
  Ag = [A zeros(n,p); -Ts*C eye(p)];
  Bg = [B; -D*Ts];
  Cg = [C zeros(p, p)];
  sysg = ss(Ag, Bg, Cg, D, Ts);
endif

[g, x, l] = lqr (sysg, q, r, s);

endfunction

%!test
%! A = [1 -1; 0 -5];
%! B = [0;1];
%! C = [1 0];
%! D = 0;
%! sys = ss(A,B,C,D);
%! [g, x, l]  = lqi(sys,eye(3),1);
%! assert(l<0,'lqi error')
