## Copyright (C) 2009-2016   Lukas F. Reichlin
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
## @deftypefn {Function File} {} nyquist (@var{sys})
## @deftypefnx {Function File} {} nyquist (@var{sys1}, @var{sys2}, @dots{}, @var{sysN})
## @deftypefnx {Function File} {} nyquist (@var{sys1}, @var{sys2}, @dots{}, @var{sysN}, @var{w})
## @deftypefnx {Function File} {} nyquist (@var{sys1}, @var{'style1'}, @dots{}, @var{sysN}, @var{'styleN'})
## @deftypefnx {Function File} {[@var{re}, @var{im}, @var{w}] =} nyquist (@var{sys})
## @deftypefnx {Function File} {[@var{re}, @var{im}, @var{w}] =} nyquist (@var{sys}, @var{w})
## Nyquist diagram of frequency response.  If no output arguments are given,
## the response is printed on the screen.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} system.  Must be a single-input and single-output (SISO) system.
## @item w
## Optional vector of frequency values.  If @var{w} is not specified,
## it is calculated by the zeros and poles of the system.
## Alternatively, the cell @code{@{wmin, wmax@}} specifies a frequency range,
## where @var{wmin} and @var{wmax} denote minimum and maximum frequencies
## in rad/s.
## @item 'style'
## Line style and color, e.g. 'r' for a solid red line or '-.k' for a dash-dotted
## black line.  See @command{help plot} for details.
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
## Version: 1.0

function [re_r, im_r, w_r] = nyquist (varargin)

  if (nargin == 0)
    print_usage ();
  endif

  [H, w, sty, leg] = __frequency_response__ ("nyquist", varargin, nargout);

  H = cellfun (@reshape, H, {[]}, {1}, "uniformoutput", false);
  re = cellfun (@real, H, "uniformoutput", false);
  im = cellfun (@imag, H, "uniformoutput", false);

  if (! nargout)
    len = numel (H);
    colororder = get (gca, "colororder");
    rc = rows (colororder);
    def_pos = arrayfun (@(k) {"-", "color", colororder(1+rem (k-1, rc), :)}, 1:len, "uniformoutput", false);
    def_neg = arrayfun (@(k) {"-.", "color", colororder(1+rem (k-1, rc), :)}, 1:len, "uniformoutput", false);
    idx = cellfun (@isempty, sty);
    sty_pos = sty_neg = sty;
    sty_pos(idx) = def_pos(idx);
    sty_neg(idx) = def_neg(idx);

    imn = cellfun (@uminus, im, "uniformoutput", false);

    pos_args = horzcat (cellfun (@horzcat, re, im, sty_pos, "uniformoutput", false){:});
    neg_args = horzcat (cellfun (@horzcat, re, imn, sty_neg, "uniformoutput", false){:});


    h = plot (pos_args{:}, neg_args{:},-1,0,'r+');
    axis ("tight")
    xlim (__axis_margin__ (xlim))
    ylim (__axis_margin__ (ylim))
    grid ("on")
    title ("Nyquist Diagram")
    xlabel ("Real Axis")
    ylabel ("Imaginary Axis")
    legend (h(1:len), leg)
  else
    re_r = re{1};
    im_r = im{1};
    w_r = w{1};
  endif

endfunction
