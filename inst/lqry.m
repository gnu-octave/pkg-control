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
## @deftypefn {Function File} {[@var{g}, @var{x}, @var{l}] =} lqry (@var{sys}, @var{q}, @var{r}, @var{s})
## Linear-quadratic regulator with output weighting.
##
## @strong{Inputs}
## @table @var
## @item sys
## Continuous or discrete-time @acronym{LTI} model (p-by-m, n states).
## @item q
## Outputs weighting matrix (p-by-p).
## @item r
## Input weighting matrix (m-by-m).
## @item s
## Optional cross term matrix (p-by-m).  If @var{s} is not specified, a zero matrix is assumed.
## @end table
##
## @strong{Outputs}
## @table @var
## @item g
## State feedback matrix (m-by-n).
## @item x
## Unique stabilizing solution of the continuous-time Riccati equation (n-by-n).
## @item l
## Closed-loop poles (n-by-1).
## @end table
##
## @strong{Equations}
## @example
## @group
## .
## x = A x + B u,   x(0) = x0
## y = C x + D u
##
##         inf
## J(x0) = INT (y' Q y  +  u' R u  +  2 y' S u)  dt
##          0
##
## L = eig (A - B*G)
## @end group
## @end example
## @seealso{lqr, dare, dlqr}
## @end deftypefn

## Author: Fabio Di Iorio <diiorio.fabio@gmail.com>
## Created: June 2024
## Version: 0.1

function [g, x, l] = lqry (sys, q, r, s = [])

  if (nargin < 3 || nargin > 4 || ~isa (sys, "lti"))
    print_usage ();
  endif

Q = q;
R = r;
N = s;

[n, m] = size(sys.b);
[p, ~] = size(sys.c);

if isempty(N)
  N=zeros(p,m);
endif

[A, B, C, D, Ts] = ssdata(sys);

[g, x, l] = lqr (sys, C'*Q*C,D'*Q*D+D'*N+N'*D+R,C'*N+C'*Q*D);

endfunction

%!test
%! A = [1 -1; 0 -5];
%! B = [0;1];
%! C = [1 0];
%! D = 0;
%! sys = ss(A,B,C,D);
%! g = lqry(sys,10,1,1);
%! assert(eig(A-B*g)<0,'lqry error')
