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
## @deftypefn {Function File} {} nichols (@var{sys})
## @deftypefnx {Function File} {} nichols (@var{sys1}, @var{sys2}, @dots{}, @var{sysN})
## @deftypefnx {Function File} {} nichols (@var{sys1}, @var{sys2}, @dots{}, @var{sysN}, @var{w})
## @deftypefnx {Function File} {} nichols (@var{sys1}, @var{'style1'}, @dots{}, @var{sysN}, @var{'styleN'})
## @deftypefnx {Function File} {[@var{mag}, @var{pha}, @var{w}] =} nichols (@var{sys})
## @deftypefnx {Function File} {[@var{mag}, @var{pha}, @var{w}] =} nichols (@var{sys}, @var{w})
## Nichols chart of frequency response.  If no output arguments are given,
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
## @seealso{bode, nyquist, sigma}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: November 2009
## Version: 1.0
##
## Author: Matthew Rhoades <mrhoa31103@gmail.com>
## Created: December 2026
## Version: 1.1
## Added Nichols Grid background
## Added Frequency Labeling on Sys Plot
##

function [mag_r, pha_r, w_r] = nichols(varargin)

  if (nargin == 0)
    print_usage ();
  endif

  [H, w, sty, sys_idx] = __frequency_response__ ("nichols", varargin, nargout);

  numsys = length (sys_idx);

  H = cellfun (@reshape, H, {[]}, {1}, "uniformoutput", false);
  mag = cellfun (@abs, H, "uniformoutput", false);
  pha = cellfun (@(H) unwrap (arg (H)) * 180 / pi, H, "uniformoutput", false);

  if (! nargout)

    ## get system names and create the legend
    leg = cell (1, numsys);
    for k = 1 : numsys
      leg{k} = inputname ( sys_idx(k) );
    endfor

    ## plot
    ## create the Nichols Grid
    ## Nichols Chart, M-Circles, and N-Circles
    ## Start with the M-Circles piece
    ## Select a figure number, clear the figure, place plot hold on so multiple
    ## circles
    ## can be drawn
    figure()
    clf
    hold on
    box on
    axis([-210 0 -24 36]);##Do not change otherwise all labels will be messed up.
    set (gca, 'xtick', -360 : 30 : 0 );
    set (gca, 'ytick', -60 : 6 : 42 );

    text(-17, -18, "-18 dB");
    text(-17, -11, "-12 dB");
    text(-17, -2, "-6 dB");
    text(-15, 3, "-4 dB");
    text(-15, 6, "-3 dB");
    text(-15, 10, "-2 dB");
    text(-15, 16, "-1 dB");
    text(-15, 22, "-.5 dB");
    text(-17, 29, "-.25 dB");
    text(-60, 32, "-.1 dB");
    text(-110, 32, "+ .1 dB");
    text(-150, 31, "+.25 dB");
    text(-138, 22.25, "+.5 dB");
    text(-140, 16.5, "+1 dB");
    text(-150, 12, "+2 dB");
    text(-172, 10.3, "+3 dB");
    text(-163, 7.7, "+4 dB");
    text(-175, 6, "+6 dB");
    text(-170, 3, "+9 dB");
    text(-175, 0,"+12 dB");
    grid

    ##Nichols chart example I have uses these gains
    dbgain = [ -18, -12, -6, -4, -3, -2, -1, -.5, -.25, -.1, 0.1, .25, .5, ...
    1, 2, 3, 4, 6, 9, 12 ];

    ##Convert to true Gains
    M = 10 .^ ( dbgain ./ 20 );
    ##create an angle vector to draw the circles
    numsegments = 1024;
    omega = linspace ( -pi, pi ,numsegments ); ## numsegments segment circle
    ##important to go from -pi to pi to prevent wrap around lines on chart.
    for i = 1 : length ( dbgain );
      ##radius of the M-circles
      radius = M(i) / ( M(i) ^2 - 1 );
      ##Generate Plot Data for M-Circle with Real Portion on Y-axis
      Xm = -M(i) ^ 2 / ( M(i) ^ 2 - 1 ) + radius * cos ( omega );
      Ym = radius * sin( omega );
      ##Create polar plot (Phase Angle vs dB Gain)
      db_plot_gain = 20 .* log10 ( sqrt (Xm .^ 2 + Ym .^ 2));
      if dbgain (i) > 0
        phase = -180 + atand ( Ym ./ Xm );
      else
        phase = atan2d ( Ym, Xm );
      endif

      for j = 1 : length ( Xm );
        if phase(j) > 0;
          phase(j) = phase(j) - 360;
        endif
      endfor
      plot ( phase, db_plot_gain, "r" );
    endfor
    ##Everything above was only to create the M Circle grid lines for the
    ##Nichols chart.

    text(-204, 21, "+2 deg", "rotation", 0);
    text(-204, 14.7, "+5 deg", "rotation", 0);
    text(-204, 9.7, "+10 deg", "rotation", 0);
    text(-204, 6, "+20 deg", "rotation", 0);
    text(-2.5, -15.842, "-2 deg", "rotation", 90);
    text(-6, -15.842, "-5 deg", "rotation", 90);
    text(-12, -15.842, "-10 deg", "rotation", 90);
    text(-23.6, -15.842, "-20 deg", "rotation", 90);
    text(-35., -15.842, "-30 deg", "rotation", 90);
    text(-68., -15.842, "-60 deg", "rotation", 90);
    text(-100., -15.842, "-90 deg", "rotation", 90);
    text(-128., -15.842, "-120 deg", "rotation", 90);
    text(-155., -15.842, "-150 deg", "rotation", 90);
    text(-182., -15.842, "-182 deg", "rotation", 90);
    text(-184., -15.842, "-185 deg", "rotation", 90);
    text(-188., -15.842, "-190 deg", "rotation", 90);
    text(-197., -15.842, "-200 deg", "rotation", 90);
    text(-205., -15.842, "-210 deg", "rotation", 90);

    ##Nichols chart example I have uses these gains
    clphase = [ -210, -178, -175, -170, -160, -150, -120, -90, -60, -20, ...
    -10, -5, -2 ];
    ##Convert to true Gains
    nphi = tand ( clphase );
    ##create an angle vector to draw the circles
    numsegments = 1024;
    omega = linspace ( -pi, pi, numsegments ); ## numsegments segment circle
    posclomega = linspace ( -pi, -pi / 2 + .035, numsegments ); ## numsegments
    ##segment circle
    ##important to go from -pi to pi to prevent wrap around lines on chart.
    for i = 1 : length ( nphi );
      ##radius of the N-circles
      radius = sqrt ( 1 / 4 + 1 / (4 * nphi(i) ^ 2 ) );
      ##Generate Plot Data for M-Circle with Real Portion on Y-axis
      if clphase(i) < 0
        Xn = -1 / 2 + radius * cos( omega );
        Yn = 1 / (2 * nphi(i)) + radius * sin ( omega );
      else
        Xn = -1 / 2 + radius * cos (posclomega);
        Yn = 1 / (2 * nphi(i)) + radius * sin ( posclomega );
      endif

      ##Create polar plot (Phase Angle vs dB Gain)
      db_plot_gain = 20 .* log10 ( sqrt( Xn .^ 2 + Yn .^ 2 ) );

      if clphase( i )>0
        phase = -180 + atand ( Yn ./ Xn );
      else
        phase = atan2d ( Yn, Xn );
      endif

      for j = 1 : length ( Xn );
        if phase(j) > 0;
          phase(j) = phase(j) - 360;
        endif
      endfor

      plot(phase, db_plot_gain, "g");

    endfor
    ##Everything above was only to create the N Circles grid lines for the Nichols chart.

    mag_db = cellfun (@mag2db, mag, "uniformoutput", false);

    plot_args = horzcat (cellfun (@horzcat, pha, mag_db, sty, "uniformoutput", false){:});
    numpoints = idivide (linspace (length ( mag_db{1} ) / 4, 3 * ...
    length ( mag_db{1} ) / 4, 10), int32 ( 1 ), "fix");
    for j = numpoints
      wstring = [num2str(w{1}(j)), " R/S"];
      text (pha{1}(j), mag_db{1}(j), wstring, "rotation", 0);
    endfor

    h11 = plot (plot_args{:});
    title ("Nichols Chart");
    xlabel ("Phase [deg]");
    ylabel ("Magnitude [dB]");
    legend (h11, leg);

  else

    ## no plotting, assign values to the output parameters
    mag_r = mag{1};
    pha_r = pha{1};
    w_r = w{1};

  endif

endfunction

%!demo
%! s = tf('s');
%! g = 1/(2*s^2+3*s+4);
%! nichols(g);
