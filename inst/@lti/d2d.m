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
## @deftypefn {Function File} {@var{sys} =} d2d (@var{sys}, @var{tsam})
## @deftypefnx {Function File} {@var{sys} =} d2d (@var{sys}, @var{tsam}, @var{method})
## @deftypefnx {Function File} {@var{sys} =} d2d (@var{sys}, @var{tsam}, @var{'prewarp'}, @var{w0})
## Resample discrete-time @acronym{LTI} model to sampling time @var{tsam}.
##
## @strong{Inputs}
## @table @var
## @item sys
## Discrete-time @acronym{LTI} model.
## @item tsam
## Desired sampling time in seconds.
## @item method
## Optional conversion method.  If not specified, default method @var{"zoh"}
## is taken.
## @table @var
## @item 'zoh'
## Zero-order hold or matrix logarithm.
## @item 'tustin', 'bilin'
## Bilinear transformation or Tustin approximation.
## @item 'prewarp'
## Bilinear transformation with pre-warping at frequency @var{w0}.
## @item 'matched'
## Matched pole/zero method.
## @end table
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## Resampled discrete-time @acronym{LTI} model with sampling time @var{tsam}.
## @end table
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2013
## Version: 0.1

function sys = d2d (sys, tsam, method = "std", w0 = 0)

  if (nargin < 2)
    print_usage ();
  endif

  tmp = d2c (sys, method, w0);

  sys = c2d (tmp, tsam, method, w0);

endfunction


%!shared num, den, z, p, k
%! H1 = zpk (0.7, 0.5, 1, 0.1);
%! H2 = d2d (H1, 0.05);
%! H3 = d2d (H2, 0.1);
%! [num, den] = tfdata (H2, "vector");
%! [z, p, k] = zpkdata (H3, "vector");
%!assert (num, [1.00000  -0.82426], 1e-4);
%!assert (den, [1.00000  -0.70711], 1e-4);
%!assert (z, 0.7, 1e-4);
%!assert (p, 0.5, 1e-4);
%!assert (k, 1.0, 1e-4);
