## Copyright (C) 2009, 2010, 2011, 2012   Lukas F. Reichlin
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
## @deftypefn {Function File} {[@var{mag}, @var{pha}, @var{w}] =} bode (@var{sys})
## @deftypefnx {Function File} {[@var{mag}, @var{pha}, @var{w}] =} bode (@var{sys}, @var{w})
## Bode diagram of frequency response.  If no output arguments are given,
## the response is printed on the screen.
##
## @strong{Inputs}
## @table @var
## @item sys
## LTI system.  Must be a single-input and single-output (SISO) system.
## @item w
## Optional vector of frequency values.  If @var{w} is not specified,
## it is calculated by the zeros and poles of the system.
## Alternatively, the cell @code{@{wmin, wmax@}} specifies a frequency range,
## where @var{wmin} and @var{wmax} denote minimum and maximum frequencies
## in rad/s.
## @end table
##
## @strong{Outputs}
## @table @var
## @item mag
## Vector of magnitude.  Has length of frequency vector @var{w}.
## @item pha
## Vector of phase.  Has length of frequency vector @var{w}.
## @item w
## Vector of frequency values used.
## @end table
##
## @seealso{nichols, nyquist, sigma}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: November 2009
## Version: 0.4

function [mag_r, pha_r, w_r] = bode2 (varargin)

  ## TODO: multiplot feature:   bode (sys1, "b", sys2, "r", ...)

  if (nargin == 0)
    print_usage ();
  endif

  [H, w] = __frequency_response_2__ (varargin, false, 0, "std", false);
  
  H = cellfun (@reshape, H, {[]}, {1}, "uniformoutput", false);
  mag = cellfun (@abs, H, "uniformoutput", false);
  pha = cellfun (@(H) unwrap (arg (H)) * 180 / pi, H, "uniformoutput", false);

  if (! nargout)
    mag_db = cellfun (@(mag) 20 * log10 (mag), mag, "uniformoutput", false);

    %w_cell = repmat ({w}, 1, numel (H));
    
    mag_args = vertcat (w, mag_db)(:);
    pha_args = vertcat (w, pha)(:);
    
    %if (isct (sys))
      xl_str = "Frequency [rad/s]";
    %else
    %  xl_str = sprintf ("Frequency [rad/s]     w_N = %g", pi / get (sys, "tsam"));
    %endif

    subplot (2, 1, 1)
    semilogx (mag_args{:})
    axis ("tight")
    ylim (__axis_margin__ (ylim))
    grid ("on")
    title (["Bode Diagram of ", inputname(1)])
    ylabel ("Magnitude [dB]")

    subplot (2, 1, 2)
    semilogx (pha_args{:})
    axis ("tight")
    ylim (__axis_margin__ (ylim))
    grid ("on")
    xlabel (xl_str)
    ylabel ("Phase [deg]")
  else
    mag_r = mag{1};
    pha_r = pha{1};
    w_r = w;
  endif

endfunction