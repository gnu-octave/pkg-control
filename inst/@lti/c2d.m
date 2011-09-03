## Copyright (C) 2009, 2011   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{sys} =} c2d (@var{sys}, @var{tsam})
## @deftypefnx {Function File} {@var{sys} =} c2d (@var{sys}, @var{tsam}, @var{method})
## Convert the continuous lti model into its discrete-time equivalent.
##
## @strong{Inputs}
## @table @var
## @item sys
## Continuous-time LTI model.
## @item tsam
## Sampling time in seconds.
## @item method
## Optional conversion method.  If not specified, default method @var{"zoh"}
## is taken.
## @table @code
## @item "zoh"
## Zero-order hold or matrix exponential.
## @item "tustin", "bilin"
## Bilinear transformation or Tustin approximation.
## @end table
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## Discrete-time LTI model.
## @end table
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.2

function sys = c2d (sys, tsam, method = "std")

  if (nargin < 2 || nargin > 3)
    print_usage ();
  endif

  if (! isa (sys, "lti"))
    error ("c2d: first argument is not an lti model");
  endif

  if (isdt (sys))
    error ("c2d: system is already discrete-time");
  endif

  if (! issample (tsam))
    error ("c2d: second argument is not a valid sample time");
  endif

  if (! ischar (method))
    error ("c2d: third argument is not a string");
  endif

  sys = __c2d__ (sys, tsam, method);
  sys.tsam = tsam;

endfunction