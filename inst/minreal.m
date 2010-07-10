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
## @deftypefn{Function File} {[@var{Amin}, @var{Bmin}, @var{Cmin}] =} ctrbf (@var{A}, @var{B}, @var{C})
## @deftypefnx{Function File} {[@var{Amin}, @var{Bmin}, @var{Cmin}] =} ctrbf (@var{A}, @var{B}, @var{C}, @var{TOL})
## Minimal realization of the system (A,B,C).
## If the system is controlable and observable, the syste es minimal.
## If the system is not controlable or/and observable, a new system created
## with the controlabe and observable subspace is equivalent 
## which means that Cco(sI-Aco)^(-1)Bco = C(sI-A)^(-1)B.
## @end deftypefn

## Author: Benjamin Fernandez <mail@benjaminfernandez.info>
## Created: 2010-07-10

function [Amin,Bmin,Cmin] = minreal(A,B,C,TOL)
  if(nargin<3 || nargin>4)
    print_usage();
  endif
  n = length(A);
  if(nargin == 3)
    TOL = n*norm(A,1)*eps;
  endif
  [Abar,Bbar,Cbar,T,K,Ac,Bc,Cc,Cnc,lc] = ctrbf(A,B,C);
  [Abar,Bbar,Cbar,T,K,Amin,Bmin,Cmin,Cno,lo] = obsvf(Ac,Bc,Cc);
  l = lc+lo;
  disp('States reduced:');
  disp(l);

endfunction
