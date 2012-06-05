## Copyright (C) 2012   Lukas F. Reichlin
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
## @deftypefn {Function File} {[@var{l}, @var{p}, @var{z}, @var{e}] =} lqe (@var{a}, @var{g}, @var{c}, @var{q}, @var{r})
## @deftypefnx {Function File} {[@var{l}, @var{p}, @var{z}, @var{e}] =} lqe (@var{a}, @var{g}, @var{c}, @var{q}, @var{r}, @var{s})
## @deftypefnx {Function File} {[@var{l}, @var{p}, @var{z}, @var{e}] =} lqe (@var{a}, @var{[]}, @var{c}, @var{q}, @var{r})
## @deftypefnx {Function File} {[@var{l}, @var{p}, @var{z}, @var{e}] =} lqe (@var{a}, @var{[]}, @var{c}, @var{q}, @var{r}, @var{s})
## Linear-quadratic estimator for discrete-time systems.
##
## @strong{Inputs}
## @table @var
## @item sys
## Continuous or discrete-time LTI model.
## @item a
## State transition matrix of continuous-time system.
## @item b
## Input matrix of continuous-time system.
## @item q
## State weighting matrix.
## @item r
## Input weighting matrix.
## @item s
## Optional cross term matrix.  If @var{s} is not specified, a zero matrix is assumed.
## @item e
## Optional descriptor matrix.  If @var{e} is not specified, an identity matrix is assumed.
## @end table
##
## @strong{Outputs}
## @table @var
## @item g
## State feedback matrix.
## @item x
## Unique stabilizing solution of the continuous-time Riccati equation.
## @item l
## Closed-loop poles.
## @end table
##
## @strong{Equations}
## @example
## @group
## .
## x = A x + B u,   x(0) = x0
##
##         inf
## J(x0) = INT (x' Q x  +  u' R u  +  2 x' S u)  dt
##          0
##
## L = eig (A - B*G)
## @end group
## @end example
## @seealso{care, dare, dlqr}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: April 2012
## Version: 0.1

function [l, p, z, e] = dlqe (a, g, c, q, r, s = [])

  if (nargin < 5 || nargin > 6)
    print_usage ();
  endif

  if (isempty (g))
    [~, p, e] = dlqr (a.', c.', q, r, s);   # dlqe (a, [], c, q, r, s), g=I
  elseif (isempty (s))
    [~, p, e] = dlqr (a.', c.', g*q*g.', r);
  else
    [~, p, e] = dlqr (a.', c.', g*q*g.', r, g*s);
  endif

  ## k computed by dlqr would be
  ## k = (r + c*p*c.') \ (c*p*a.' + s.')
  ## such that  l = a \ k.'
  ## what about the s term?
    
  l = p*c.' / (c*p*c.' + r);
  ## z = p - p*c.' / (c*p*c.' + r) * c*p;
  z = p - l*c*p;

endfunction
