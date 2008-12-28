## Copyright (C) 1996, 2000, 2003, 2004, 2005, 2007
##               Auburn University. All rights reserved.
##
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {} dgram (@var{a}, @var{b})
## Return controllability gramian of discrete time system
## @iftex
## @tex
## $$ x_{k+1} = ax_k + bu_k $$
## @end tex
## @end iftex
## @ifinfo
## @example
##   x(k+1) = a x(k) + b u(k)
## @end example
## @end ifinfo
## 
## @strong{Inputs}
## @table @var
## @item a
## @var{n} by @var{n} matrix
## @item b
## @var{n} by @var{m} matrix
## @end table
##
## @strong{Output}
## @table @var
## @item m 
## @var{n} by @var{n} matrix, satisfies
## @iftex
## @tex
## $$ ama^T - m + bb^T = 0 $$
## @end tex
## @end iftex
## @ifinfo
## @example
##  a m a' - m + b*b' = 0
## @end example
## @end ifinfo
## @end table
## @end deftypefn

## Author: A. S. Hodel <a.s.hodel@eng.auburn.edu>
## Created: July 1995

function m = dgram (a, b)

  if (nargin != 2)
    print_usage ();
  else

    ## this is not efficient, but it is simple to mantain
    ##
    ## efficiency is not important because this function is here only
    ## for backward compatibility
    c = zeros (1, length (a)); ## it doesn't matter what the value of c is
    d = zeros (rows (c), columns (b)); ## it doesn't matter what the value of d is
    Ts = 0.1; ## Ts != 0
    sys = ss (a, b, c, d, Ts);
    m = gram (sys, 'c');
  endif

endfunction


%!test
%! a = [-1 0 0; 1/2 1 0; 1/2 0 -1] / 2;
%! b = [1 0; 0 -1; 0 1];
%! Wc = dgram (a, b);
%! assert (a * Wc * a' - Wc + b * b', zeros (size (a)), 1e-12)