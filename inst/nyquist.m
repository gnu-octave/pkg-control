## Copyright (C) 2009 - 2010   Lukas F. Reichlin
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
## @deftypefn {Function File} {[@var{re}, @var{im}, @var{w}] =} nyquist (@var{sys})
## @deftypefnx {Function File} {[@var{re}, @var{im}, @var{w}] =} nyquist (@var{sys}, @var{w})
## Nyquist diagram of frequency response. If no output arguments are given,
## the response is printed on the screen.
##
## @strong{Inputs}
## @table @var
## @item sys
## LTI system. Must be a single-input and single-output (SISO) system.
## @item w
## Optional vector of frequency values. If @var{w} is not specified, it
## is calculated by the zeros and poles of the system.
## @end table
##
## @strong{Outputs}
## @table @var
## @item re
## Vector of real parts. Has length of frequency vector @var{w}.
## @item im
## Vector of imaginary parts. Has length of frequency vector @var{w}.
## @item w
## Vector of frequency values used.
## @end table
##
## @seealso{bode, nichols, sigma}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: November 2009
## Version: 0.2

function [re_r, im_r, w_r] = nyquist (sys, w = [])

  ## TODO: multiplot feature:   nyquist (sys1, "b", sys2, "r", ...)

  if (nargin == 0 || nargin > 2)
    print_usage ();
  endif

  [H, w] = __frequency_response__ (sys, w, false, 0, "ext");

  H = reshape (H, [], 1);
  re = real (H);
  im = imag (H);

  if (! nargout)
    ax_vec = __axis2dlim__ ([re, im; re, -im]);

    plot (re, im, "b", re, -im, "r")
    axis (ax_vec)
    grid ("on")
    title ("Nyquist Diagram")
    xlabel ("Real Axis")
    ylabel ("Imaginary Axis")
  else
    re_r = re;
    im_r = im;
    w_r = w;
  endif

endfunction