## Copyright (C) 2010 Benjamin Fernandez
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn{Function File} {[@var{Abar}, @var{Bbar}, @var{Cbar}, @var{T}, @var{K}] =} obsvf (@var{A}, @var{B}, @var{C})
## @deftypefnx{Function File} {[@var{Abar}, @var{Bbar}, @var{Cbar}, @var{T}, @var{K}] =} obsvf (@var{A}, @var{B}, @var{C}, @var{TOL})
## If Ob=obsv(A,C) has rank r <= n = SIZE(A,1), then there is a 
## similarity transformation Tc such that To = [t1;t2] where t1 is c
## and t2 is orthogonal to t1
##
## @example
## @group
## Abar = To \ A * To ,  Bbar = To \ B ,  Cbar = C * To
## @end group
## @end example
## 
## and the transformed system has the form
##
## @example
## @group
##        | Ao     0 |           | Bo  |
## Abar = |----------|,   Bbar = | --- |,  Cbar = [Co | 0 ].
##        | A21   Ano|           | Bno |
## @end group
## @end example
##                                                      
## where (Ao,Bo) is observable, and Co(sI-Ao)^(-1)Bo = C(sI-A)^(-1)B. And 
## system is detectable if Ano has no eigenvalues in the right 
## half plane. The last output K is a vector of length n containing the 
## number of observable states.
## @end deftypefn

## Author: Benjamin Fernandez <benjas@benjas-laptop>
## Created: 2010-05-02

function [Abar, Bbar, Cbar, T, K] = obsvf (A, B, C, TOL)

  if (nargin<3 || nargin>4)
    print_usage ();
  end

  n = length (A);

  if (nargin==3)
    TOL = n*norm (A, 1)*eps;
  end

  Ob         = obsv (A, C);
  [nro, nco] = size (Ob);
  rob        = rank (Ob);
  lr         = nco-rob;
  [U, S, V]  = svd (Ob);
  K          = V(:, 1:rob);  # Basis raw space
  T          = V;  # [c;orth(c)];
  Abar       = T\A*T;
  Bbar       = T\B;
  Cbar       = C*T;
  Ano        = Abar(n-(n-rob)+1:n, n-(n-rob)+1:n)
  Ao         = Abar(1:rob, 1:rob);
  Bo         = Bbar(1:rob, :);
  Co         = Cbar(:, 1:rob);

endfunction