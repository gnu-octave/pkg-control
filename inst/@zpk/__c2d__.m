## Copyright (C) 2026  Mitchell Thompkins <mitchell.thompkins@pm.me>
##
## This file is part of the control package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not,
## see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## Convert the continuous ZPK model into its discrete-time equivalent.
## The matched method maps each stored pole and zero directly via exp (s*tsam),
## avoiding the ill-conditioned polynomial round-trip of the TF representation.

## Author: Mitchell Thompkins <mitchell.thompkins@pm.me>
## Created: June 2026
## Version: 0.1

function sys = __c2d__ (sys, tsam, method = "zoh", w0 = 0)

  if (strncmpi (method, "m", 1))    # "matched"

    if (! issiso (sys))
      error ("zpk: c2d: require SISO system for matched pole/zero method");
    endif

    z_c = sys.z{1};
    p_c = sys.p{1};
    k_c = sys.k;

    p_d = exp (p_c * tsam);
    z_d = exp (z_c * tsam);

    if (any (! isfinite (p_d)) || any (! isfinite (z_d)))
      error ("zpk: c2d: discrete-time poles and zeros are not finite");
    endif

    ## continuous-time zeros at infinity are mapped to -1 in discrete-time
    ## except for one.  for non-proper transfer functions, no zeros at -1 are added.
    np = length (p_c);              # number of poles
    nz = length (z_c);              # number of finite zeros, np-nz number of infinite zeros
    z_d = vertcat (z_d, repmat (-1, np-nz-1, 1));

    ## the discrete-time gain k_d is matched at frequency w_c to continuous-time
    ## gain k_c.  dc gain is taken (w_c=0) unless there are continuous-time
    ## poles/zeros near the imaginary axis at j*w_c.  gain is evaluated on the
    ## imaginary axis (s=j*w_c) and unit circle (z=exp(j*w_c*tsam)) so that
    ## |H_d(exp(j*w_c*tsam))| = |H_c(j*w_c)| holds in the frequency domain.
    w_c = 0;                        # start at dc
    tol = sqrt (eps);               # poles/zeros within tol of j*w_c are avoided
    while (any (abs ([p_c; z_c] - 1j*w_c) < tol))
      w_c += 0.1 / tsam;
    endwhile
    w_d = exp (1j * w_c * tsam);
    k_d = real (k_c * prod (1j*w_c - z_c) / prod (1j*w_c - p_c) * prod (w_d - p_d) / prod (w_d - z_d));

    sys.z{1} = z_d;
    sys.p{1} = p_d;
    sys.k = k_d;

  else
    ## zoh/foh/tustin/prewarp/impulse are not per-root maps. their natural
    ## representation is state-space.  convert back to zpk for type consistency.
    sys = zpk (__c2d__ (ss (sys), tsam, method, w0));
  endif

endfunction


%!test
%! ## single pole: p_d = exp(p_c * Ts)
%! sys = zpk ([], [-1], 1);
%! sys_d = c2d (sys, 0.1, 'matched');
%! assert (isa (sys_d, 'zpk'));
%! [~, p_d] = zpkdata (sys_d, 'v');
%! assert (p_d, exp (-0.1), 1e-14);

%!test
%! ## 25 poles clustered near the imaginary axis: matched c2d keeps them
%! ## stable and maps each one exactly to exp(p*Ts) with no rounding error
%! N = 25; Ts = 1/1000;
%! p_s = (-0.001 + 1i * linspace (1, 25, N))' * 2*pi*4;
%! sys_d = c2d (zpk ([], p_s, 1), Ts, 'matched');
%! [~, p_d] = zpkdata (sys_d, 'v');
%! assert (all (abs (p_d) < 1));
%! assert (max (abs (p_d - exp (p_s * Ts))), 0, 1e-12);

%!test
%! ## non-matched methods still work on zpk and return zpk
%! sys_d = c2d (zpk ([], [-1], 1), 0.1, 'zoh');
%! assert (isa (sys_d, 'zpk'));
%! assert (get (sys_d, 'tsam'), 0.1);
