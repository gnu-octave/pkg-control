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
## @deftypefn {Function File} {@var{k} =} dcgain (@var{sys})
## Compute the DC gain of @acronym{LTI} system.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} system created by tf(), ss(), dss(), etc.
## @end table
##
## @strong{Outputs}
## @table @var
## @item k
## DC gain matrix. For a system with m inputs and p outputs, the array @var{k}
## has dimensions [p, m].
## @end table
##
## Transfer function for a continuous state space system (A,B,C,D)
## G(s) = C * inv(s*I - A) * B + D
##
## DC Gain: evaluate G(s) as s -> 0:
## k = C * inv(-A) * B + D
##
## Transfer function for a discrete state space system (A,B,C,D,T)
## G(z) = C * inv(z*I - A) * B + D
##
## DC Gain: evaluate G(z) as z -> 1:
## k = C * inv(I-A) * B + D
##
## @strong{Example}
## @example
## G = Transfer function 'G' from input 'u1' to output ...
##
##                1          
## y1:  ---------------------
##      s^3 + 2 s^2 + 3 s + 4
##
##
## octave:1> K = dcgain(G)
##
## K =  0.25000
## @end example

## @seealso{freqresp,tf,ss,dss}
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

%!assert( dcgain( tf(1,[1,1]) )                   , 1   )
%!assert( dcgain( tf(2,[1,1]) )                   , 2   )
%!assert( dcgain( ss([0,1;-2,-3],[0;1],[1,0],0) ) , 0.5 )
