## Copyright (C) 2010 Benjamin Fernandez <mail@benjaminfernandez.info>
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
## @deftypefn{Function File} {[@var{Abar}, @var{Bbar}, @var{Cbar}, @var{T}, @var{K}] =} ctrbf (@var{A}, @var{B}, @var{C})
## @deftypefnx{Function File} {[@var{Abar}, @var{Bbar}, @var{Cbar}, @var{T}, @var{K}] =} ctrbf (@var{A}, @var{B}, @var{C}, @var{TOL})
## If Co=ctrb(A,B) has rank r <= n = SIZE(A,1), then there is a 
## similarity transformation Tc such that Tc = [t1 t2] where t1
## is the controllable subspace and t2 is orthogonal to t1
##
## @example
## @group
## Abar = Tc \ A * Tc ,  Bbar = Tc \ B ,  Cbar = C * Tc
## @end group
## @end example
##
## and the transformed system has the form
##
## @example
## @group
##        | Ac    A12|           | Bc |
## Abar = |----------|,   Bbar = | ---|,  Cbar = [Cc | Cnc].
##        | 0     Anc|           |  0 |
## @end group
## @end example
##                                     
## where (Ac,Bc) is controllable, and Cc(sI-Ac)^(-1)Bc = C(sI-A)^(-1)B.
## and the system is stabilizable if Anc has no eigenvalues in
## the right half plane. The last output K is a vector of length n
## containing the number of controllable states.
## @end deftypefn

## Author: Benjamin Fernandez <benjas@benjas-laptop>
## Created: 2010-04-30

function [Abar, Bbar, Cbar, T, K] = ctrbf (A, B, C, TOL)

  if (nargin < 3 || nargin > 4)
    print_usage ();
  endif

  if (nargin == 3)
    TOL = length (A) * norm (A,1) * eps;
  endif

  Co         = ctrb (A, B);
  [nrc, ncc] = size (Co);
  rco        = rank (Co, TOL);
  lr         = nrc - rco;
  [U, S, V]  = svd (Co);
  K          = U(:, 1:rco);  # Basis column space
  T          = U;            # [B orth(B)]
  Abar       = T \ A * T;
  Bbar       = T \ B;
  Cbar       = C * T;
  
endfunction