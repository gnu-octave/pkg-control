## Copyright (C) 2011   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{sys} =} d2c (@var{sys})
## @deftypefnx {Function File} {@var{sys} =} d2c (@var{sys}, @var{method})
## Convert the discrete lti model into its continuous-time equivalent.
##
## @strong{Inputs}
## @table @var
## @item sys
## Discrete-time LTI model.
## @item method
## Optional conversion method.  If not specified, default method @var{"zoh"}
## is taken.
## @table @code
## @item "zoh"
## Zero-order hold or matrix logarithm.
## @item "tustin", "bilin"
## Bilinear transformation or Tustin approximation.
## @end table
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## Continuous-time LTI model.
## @end table
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2011
## Version: 0.1

function sys = d2c (sys, method = "std")

  if (nargin == 0 || nargin > 2)
    print_usage ();
  endif

  if (! isa (sys, "lti"))
    error ("d2c: first argument is not an lti model");
  endif

  if (isct (sys))
    error ("d2c: system is already continuous-time");
  endif

  if (! ischar (method))
    error ("c2d: second argument is not a string");
  endif

  sys = __d2c__ (sys, sys.tsam, method);
  sys.tsam = 0;

endfunction