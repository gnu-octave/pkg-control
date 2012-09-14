## Copyright (C) 2009, 2010, 2012   Lukas F. Reichlin
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
## @deftypefn {Function File} {[@var{re}, @var{im}, @var{w}] =} nyquist (@var{sys})
## @deftypefnx {Function File} {[@var{re}, @var{im}, @var{w}] =} nyquist (@var{sys}, @var{w})
## Nyquist diagram of frequency response.  If no output arguments are given,
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
## @item re
## Vector of real parts.  Has length of frequency vector @var{w}.
## @item im
## Vector of imaginary parts.  Has length of frequency vector @var{w}.
## @item w
## Vector of frequency values used.
## @end table
##
## @seealso{bode, nichols, sigma}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: November 2009
## Version: 0.4

function [re_r, im_r, w_r] = nyquist2 (varargin)

  if (nargin == 0)
    print_usage ();
  endif

  [H, w] = __frequency_response_2__ (varargin, false, 0, "ext");

  H = cellfun (@reshape, H, {[]}, {1}, "uniformoutput", false);
  re = cellfun (@real, H, "uniformoutput", false);
  im = cellfun (@imag, H, "uniformoutput", false);

  if (! nargout)
    tmp = cellfun (@isa, varargin, {"lti"});
    sys_idx = find (tmp);
    tmp = cellfun (@ischar, varargin);
    style_idx = find (tmp);

    len = numel (H);  
    plot_args_pos = {};
    plot_args_neg = {};
    legend_args = cell (len, 1);
    colororder = get (gca, "colororder");
    rc = rows (colororder);

    for k = 1:len
      col = colororder(1+rem (k-1, rc), :);
      plot_args_pos = cat (2, plot_args_pos, re{k}, im{k}, {"-", "color", col});
      plot_args_neg = cat (2, plot_args_neg, re{k}, -im{k}, {"-.", "color", col});
      legend_args{k} = inputname(sys_idx(k));
    endfor

    h = plot (plot_args_pos{:}, plot_args_neg{:});
    axis ("tight")
    xlim (__axis_margin__ (xlim))
    ylim (__axis_margin__ (ylim))
    grid ("on")
    title ("Nyquist Diagram")
    xlabel ("Real Axis")
    ylabel ("Imaginary Axis")
    legend (h(1:len), legend_args)
    % legend (h(1:2:2*len), legend_args)
  else
    re_r = re{1};
    im_r = im{1};
    w_r = w{1};
  endif

endfunction