## Copyright (C) 2013-2015   Thomas Vasileiou
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
## @deftypefn {Function File} {@var{sys} =} thiran (@var{tau}, @var{tsam})
## Approximation of continuous-time delay using a discrete-time
## allpass Thiran filter. 
##
## Thiran filters can approximate continuous-time delays that are
## non-integer multiples of the sampling time (fractional delays).
## This approximation gives a better matching of the phase shift
## between the continuous- and the discrete-time system.
## If there is no fractional part in the delay, then the standard
## discrete-time delay representation is used. 
##
## @strong{Inputs}
## @table @var
## @item tau
## A continuous-time delay, given in time units (seconds).
## @item tsam
## The sampling time of the resulting Thiran filter.
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## Transfer function model of the resulting filter.
## The order of the filter is determined automatically.
## @end table
##
## @strong{Example}
## @example
## @group
## octave:1> sys = thiran (1.33, 0.5)
##
## Transfer function 'sys' from input 'u1' to output ...
## 
##       0.003859 z^3 - 0.03947 z^2 + 0.2787 z + 1
##  y1:  -----------------------------------------
##        z^3 + 0.2787 z^2 - 0.03947 z + 0.003859 
## 
## Sampling time: 0.5 s
## Discrete-time model.
## @end group
## @end example
## @example
## @group
## octave:2> sys = thiran (1, 0.5)
##
## Transfer function 'sys' from input 'u1' to output ...
## 
##        1
##  y1:  ---
##       z^2 
## 
## Sampling time: 0.5 s
## Discrete-time model.
## @end group
## @end example
##
## @seealso{absorbdelay, pade}
## @end deftypefn

## Author: Thomas Vasileiou <thomas-v@wildmail.com>
## Created: January 2013
## Version: 0.1

function sys = thiran (del, tsam)

  ## check args
  if (nargin != 2)
    print_usage ();
  endif

  if (! issample (del, 0))
    error ("thiran: the delay parameter 'tau' must be a non-negative scalar.");
  endif
  
  if (! issample (tsam))
    error ("thiran: the second parameter 'tsam' is not a valid sampling time.");
  endif
  
  if (del == 0)
    sys = tf (1);
    return;
  endif
  
  ## find fractional and discrete delay
  N = floor (del/tsam + eps);       # put eps or sometimes it misses
  d = del - N*tsam;
  
  ## check if we do need a thiran filter
  if (d/tsam < eps)
    sys = tf (1, [1, zeros(1, N)], tsam);
  else
    ## make filter order ~ del to ensure stability
    N = N + 1;                      # order of the filter
    d = del/tsam;
    tmp = ((d-N+(0:N).') * ones (1,N)) ./ (d-N + bsxfun (@plus, 1:N, (0:N).'));
    a = horzcat (1, (-1).^(1:N) .* bincoeff (N, 1:N) .* prod (tmp));
    sys = tf (fliplr (a), a, tsam);
  endif

endfunction


%!shared num, den, expc
%! expc = [1, 0.5294, -0.04813, 0.004159];
%! sys = thiran (2.4, 1);
%! [num, den] = tfdata (sys, "vector");
%!assert (den, expc, 1e-4);
%!assert (num, fliplr (expc), 1e-4);

%!shared num, den
%! sys = thiran (0.5, 0.1);
%! [num, den] = tfdata (sys, "vector");
%!assert (den, [1, 0, 0, 0, 0, 0], eps);
%!assert (num, 1, eps);

%!error (thiran (-1, 1));
%!error (thiran (1, -1));
%!error (thiran ([1 2 3], 1));
