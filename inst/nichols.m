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
## @deftypefn {Function File} {[@var{mag}, @var{pha}, @var{w}] =} nichols (@var{sys})
## @deftypefnx {Function File} {[@var{mag}, @var{pha}, @var{w}] =} nichols (@var{sys}, @var{w})
## Nichols chart of LTI model's frequency response.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: November 2009
## Version: 0.1

function [mag_r, pha_r, w_r] = nichols (sys, w = [])

  if (nargin == 0 || nargin > 2)
    print_usage ();
  endif

  [H, w] = __getfreqresp__ (sys, w, false, 0, "ext");

  H = H(:);
  mag = abs (H);
  pha = unwrap (arg (H)) * 180 / pi;

  if (! nargout)
    mag_db = 20 * log10 (mag);
    ax_vec = __axis2dlim__ ([pha(:), mag_db(:)]);
    
    plot (pha, mag_db)
    axis (ax_vec)
    grid ("on")
    title ("Nichols Chart")
    xlabel ("Phase [deg]")
    ylabel ("Magnitude [dB]")
  else
    mag_r = mag;
    pha_r = pha;
    w_r = w;
  endif

endfunction