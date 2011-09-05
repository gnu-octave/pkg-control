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
## Convert the continuous SS model into its discrete-time equivalent.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.2

function sys = __c2d__ (sys, tsam, method = "zoh")

  switch (method)
    case {"zoh", "std"}
      [sys.a, sys.b, sys.c, sys.d, sys.e] = __dss2ss__ (sys.a, sys.b, sys.c, sys.d, sys.e);
      [n, m] = size (sys.b);       # n: states, m: inputs
      ## TODO: use SLICOT MB05OD
      tmp = expm ([sys.a*tsam, sys.b*tsam; zeros(m, n+m)]);
      sys.a = tmp (1:n, 1:n);      # F
      sys.b = tmp (1:n, n+(1:m));  # G

    case {"tustin", "bilin"}
      [sys.a, sys.b, sys.c, sys.d, sys.e] = __dss2ss__ (sys.a, sys.b, sys.c, sys.d, sys.e);
      [sys.a, sys.b, sys.c, sys.d] = slab04md (sys.a, sys.b, sys.c, sys.d, 1, 2/tsam, false);
      ## TODO: descriptor case

    ## TODO: case "prewarp"

    otherwise
      error ("ss: c2d: %s is an invalid or missing method", method);
  endswitch

endfunction
