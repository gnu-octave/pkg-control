## Copyright (C) 2009-2016   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{K} =} dcgain (@var{sys})
## Compute the DC gain of @acronym{LTI} systems given as transfer function.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} system created by tf(), ss(), dss(), etc.
## @end table
##
## @strong{Outputs}
## @table @var
## @item K
## DC gain matrix. For a system with m inputs and p outputs, the array @var{k}
## has dimensions [p, m].
## @end table
##
## @xref{@@tf/dcgain} for the method used for systems represented by a transfer function.
##
## For a continuous-time state space system @math{(A,B,C,D)}, the gain is given by
## @tex
## $$ K = G(s=0) = \left. C \, (s\,I - A)^{-1} B + D\,\right|_{s=0}
##      = C \, A^{-1} B + D $$
## @end tex
## @ifnottex
## @example
## @group
##   K = G(s=0) = C * inv(A) * B + D
## @end group
## @end example
## @end ifnottex
##
## In the discrete-time case the gain of a state space system @math{(A,B,C,D,T)}
## is given by
## @tex
## $$ K = G(z=1) = \left. C \, (z\,I - A)^{-1} B + D\,\right|_{z=1}
##      = C \, (I-A)^{-1} B + D $$
## @end tex
## @ifnottex
## @example
## @group
##   K = G(z=1) = C * inv(I - A) * B + D
## @end group
## @end example
## @end ifnottex
##
## @strong{Example}
## @example
## @group
## A = [0,1,0;0,0,1;-1,-3,-3];
## B = [0 0;0,1;1,1];
## C = [1,0 0;0,1,1];
## D = [2,0;0,0];
## K = dcgain (ss (A,B,C,D))
##
## K =
##   3   4
##   0  -1
## @end group
## @end example
##
## @seealso{@@tf/dcgain,@@lti/freqresp,@@tf/tf,@@ss/ss,dss}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Author: Geraint Paul Bevan <geraint.bevan@gcu.ac.uk>
## Created: October 2009
## Version: 0.1

function gain = dcgain (sys)

  if (nargin != 1)  # sys is always an LTI model
    print_usage ();
  endif

  gain = __freqresp__ (sys, 0);

endfunction

%!assert( dcgain( ss([0,1;-2,-3],[0;1],[1,0],0) ) , 0.5 )
%!assert( dcgain( ss([0,1,0;0,0,1;-1,-3,-3],[0 0;0,1;1,1],[1,0 0;0,1,1],[2,0;0,0]) ) , [3,4;0,-1])
