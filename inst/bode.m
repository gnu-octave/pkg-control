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
## @deftypefn {Function File} {} bode (@var{sys})
## @deftypefnx {Function File} {} bode (@var{sys1}, @var{sys2}, @dots{}, @var{sysN})
## @deftypefnx {Function File} {} bode (@var{sys1}, @var{sys2}, @dots{}, @var{sysN}, @var{w})
## @deftypefnx {Function File} {} bode (@var{sys1}, @var{'style1'}, @dots{}, @var{sysN}, @var{'styleN'})
## @deftypefnx {Function File} {[@var{mag}, @var{pha}, @var{w}] =} bode (@var{sys})
## @deftypefnx {Function File} {[@var{mag}, @var{pha}, @var{w}] =} bode (@var{sys}, @var{w})
## Bode diagram of frequency response.  If no output arguments are given,
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
## Version: 1.0

function [mag_r, pha_r, w_r] = bode (varargin)

  if (nargin == 0)
    print_usage ();
  endif

  [H, w, sty, sys_idx, H_auto, w_auto] = __frequency_response__ ("bode", varargin, nargout);

  H = cellfun (@reshape, H, {[]}, {1}, "uniformoutput", false);
  H_auto = cellfun (@reshape, H_auto, {[]}, {1}, "uniformoutput", false);
  mag = cellfun (@abs, H, "uniformoutput", false);
  pha = cellfun (@(H) unwrap (arg (H)) * 180 / pi, H, "uniformoutput", false);
  pha_auto = cellfun (@(H_auto) unwrap (arg (H_auto)) * 180 / pi, H_auto, "uniformoutput", false);
  numsys = length (sys_idx);

  ## check for poles and zeroes at the origin for each of the numsys systems
  ## and compare this asymptotic value to initial phase for the full auto
  ## w-range, which is supposed to start at sufficiently small values fo w
  initial_phase = cellfun(@(x) x(1), pha_auto);
  for h = 1:numsys
    sys1 = varargin{sys_idx(h)};
    [num,den] = tfdata (sys1,'v');
    n_poles_at_origin = sum (roots (den) == 0);
    n_zeros_at_origin = sum (roots (num) == 0);
    asymptotic_low_freq_phase = (n_zeros_at_origin - n_poles_at_origin)*90;
    pha_auto(h)={cell2mat(pha_auto(h)) + round((asymptotic_low_freq_phase - initial_phase(h))/360)*360};
  endfor

  ## check if possibly given w-range is at higher frequencies and provide
  ## missing "history" for the the phase
  initial_phase = cellfun(@(x) x(1), pha);
  for h = 1:numsys
    if (length (w{h}) != length (w_auto{h})) || any (w{h} - w_auto{h})
      ## the w-range to use is not the auto (full) range, thus, make sure
      ## the inforamtion on the phase form beginning  at small frequencies
      ## is used to determine the correct pahse at the beginning of the desired
      ## w range
      if w{h}(1) > w_auto{h}(1)
        ## w range starts at higher frequencies; get the nearest w in the auto
        ## range and add +/-360 degrees minimizing the differenze betwenn the
        ## phases
        pha_idx = length (find (w_auto{h} < w{h}(1)));
        pha_cmp = pha_auto{h}(pha_idx);
        pha(h)={cell2mat(pha(h)) + round((pha_cmp - initial_phase(h))/360)*360};
      endif
    else
      ## w range is identical to auto range, just use the phase for the
      ## auto range
      pha(h) = {pha_auto(h)};
    endif
  endfor


  if (! nargout)

    ## get system names and create the legend
    leg = cell (1, numsys);
    for k = 1:numsys
      leg{k} = inputname (sys_idx(k));
      if (isempty (leg{k}))
        leg{k} = sprintf("Sys %d", k);
      endif
    endfor

    ## plot
    mag_db = cellfun (@mag2db, mag, "uniformoutput", false);

    mag_args = horzcat (cellfun (@horzcat, w, mag_db, sty, "uniformoutput", false){:});
    pha_args = horzcat (cellfun (@horzcat, w, pha, sty, "uniformoutput", false){:});

    hax1 = subplot (2, 1, 1)
    semilogx (mag_args{:})
    axis ("tight")
    ylim (__axis_margin__ (ylim))
    grid ("on")
    title ("Bode Diagram")
    ylabel ("Magnitude [dB]")

    subplot (2, 1, 2)
    semilogx (pha_args{:})
    axis ("tight")
    ylim (__axis_margin__ (ylim))
    grid ("on")
    xlabel ("Frequency [rad/s]")
    ylabel ("Phase [deg]")

    ## Make top axes current for ML compatibility
    set (gcf, "currentaxes", hax1);

  else

    ## no plotting, assign values to the output parameters
    mag_r = mag{1};
    pha_r = pha{1};
    if (iscell (pha_r)) # fix possible cell in cell
      pha_r = pha_r{1};
    endif
    w_r = w{1};

  endif

endfunction

%!demo
%! s = tf('s');
%! g = 1/(2*s^2+3*s+4);
%! bode(g);

%!test
%! s = tf('s');
%! K = 1;
%! T = 2;
%! g = K/(1+T*s);
%! [mag phas w] = bode(g);
%! mag_dB = 20*log10(mag);
%! index = find(mag_dB < -3,1);
%! w_cutoff = w(index);
%! assert (1/T, w_cutoff, eps);
