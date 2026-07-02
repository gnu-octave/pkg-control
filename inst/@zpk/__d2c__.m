## Copyright (C) 2026  Mitchell Thompkins
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
## Convert the discrete ZPK model into its continuous-time equivalent.

## Author: Mitchell Thompkins <mitchell.thompkins@pm.me>
## Created: June 2026
## Version: 0.1

function sys = __d2c__ (sys, tsam, method = "zoh", w0 = 0)

  if (strncmpi (method, "m", 1))    # "matched"

    if (! issiso (sys))
      error ("zpk: d2c: require SISO system for matched pole/zero method");
    endif

    z_d = sys.z{1};
    p_d = sys.p{1};
    k_d = sys.k;

    if (any (abs (p_d) < eps) || any (abs (z_d) < eps))
      error ("zpk: d2c: discrete-time poles and zeros at 0 not supported because log(0) is -Inf");
    endif

    z_d_orig = z_d;
    z_d(abs (z_d+1) < sqrt (eps)) = [];    # remove zeros added at -1 by c2d

    p_c = log (p_d) / tsam;
    z_c = log (z_d) / tsam;

    w_c = 0;
    w_d = 1;
    tol = sqrt (eps);
    while (any (abs ([p_d; z_d_orig] - w_d) < tol))
      w_c += 0.1 / tsam;
      w_d = exp (w_c * tsam);
    endwhile
    k_c = real (k_d * prod (w_d - z_d_orig) / prod (w_d - p_d) * prod (w_c - p_c) / prod (w_c - z_c));

    sys.z{1} = z_c;
    sys.p{1} = p_c;
    sys.k = k_c;

  else
    sys = zpk (__d2c__ (ss (sys), tsam, method, w0));
  endif

endfunction


%!test
%! ## matched round trip is the identity on poles/zeros
%! sys = zpk ([-3], [-1; -2], 4);
%! sysd = c2d (sys, 0.1, 'matched');
%! sysc = d2c (sysd, 'matched');
%! [z, p, k] = zpkdata (sysc, 'v');
%! assert (sort (real (p)), [-2; -1], 1e-10);
%! assert (real (z), -3, 1e-10);
%! assert (k, 4, 1e-10);
