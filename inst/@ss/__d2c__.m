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
## Convert the discrete SS model into its continuous-time equivalent.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2011
## Version: 0.1

function sys = __d2c__ (sys, tsam, method = "zoh")

  switch (method)
    case {"zoh", "std"}
      if (! isempty (sys.e))
        if (rcond (sys.e) < eps)
          error ("ss: d2c: zero-order hold method requires proper system");
        else
          sys.a = sys.e \ sys.a;
          sys.b = sys.e \ sys.b;
          sys.e = [];              # require ordinary state-space model
        endif
      endif

      error ("ss: d2c: zoh method not implemented yet");
      [n, m] = size (sys.b);       # n: states, m: inputs

      ## TODO: use SLICOT MB05OD
      tmp = expm ([sys.a*tsam, sys.b*tsam; zeros(m, n+m)]);

      sys.a = tmp (1:n, 1:n);      # F
      sys.b = tmp (1:n, n+(1:m));  # G

    case {"tustin", "bilin"}
      if (! isempty (sys.e))
        if (rcond (sys.e) < eps)
          error ("ss: d2c: tustin method requires proper system");
        else
          sys.a = sys.e \ sys.a;
          sys.b = sys.e \ sys.b;
          sys.e = [];              # require ordinary state-space model
        endif
      endif

      [sys.a, sys.b, sys.c, sys.d] = slab04md (sys.a, sys.b, sys.c, sys.d, 1, 2/tsam, true);

    otherwise
      error ("ss: d2c: %s is an invalid or missing method", method);

  endswitch

endfunction