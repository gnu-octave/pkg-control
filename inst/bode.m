## Copyright (C) 2009   Lukas F. Reichlin
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
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{mag}, @var{pha}, @var{w}] =} bode (@var{sys})
## @deftypefnx {Function File} {[@var{mag}, @var{pha}, @var{w}] =} bode (@var{sys}, @var{w})
## Bode diagram of LTI model's frequency response.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: November 2009
## Version: 0.1

function [mag_r, pha_r, w_r] = bode (sys, w = [])

  if (nargin == 0 || nargin > 2)
    print_usage ();
  endif

  [H, w] = __getfreqresp__ (sys, w, false, 0, "std");

  H = H(:);
  mag = abs (H);
  pha = unwrap (arg (H)) * 180 / pi;

  if (! nargout)
    mag_db = 20 * log10 (mag);

    wv = [min(w), max(w)];
    ax_vec_mag = __axis2dlim__ ([w(:), mag_db(:)]);
    ax_vec_mag(1:2) = wv;
    ax_vec_pha = __axis2dlim__ ([w(:), pha(:)]);
    ax_vec_pha(1:2) = wv;

    if (isct (sys))
      xl_str = "Frequency [rad/s]";
    else
      xl_str = sprintf ("Frequency [rad/s]     w_N = %g", pi / get (sys, "tsam"));
    endif

    subplot (2, 1, 1)
    semilogx (w, mag_db)
    axis (ax_vec_mag)
    grid ("on")
    title ("Bode Diagram")
    ylabel ("Magnitude [dB]")

    subplot (2, 1, 2)
    semilogx (w, pha)
    axis (ax_vec_pha)
    grid ("on")
    xlabel (xl_str)
    ylabel ("Phase [deg]")
  else
    mag_r = mag;
    pha_r = pha;
    w_r = w;
  endif

endfunction