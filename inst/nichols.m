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
## Nichols chart of frequency response, plots gain over phase.
## If no output arguments are given, the response is printed on
## the screen together with a grid of contours of constant close-loop
## magnitude and phase.
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
    for k = 1:numsys
      leg{k} = inputname (sys_idx(k));
      if (isempty (leg{k}))
        leg{k} = sprintf("Sys %d", k);
      endif
    endfor

    ## plot
    ## create the Nichols Grid
    ## Nichols Chart, M-Circles, and N-Circles
    ## Start with the M-Circles piece
    ## Select a figure number, clear the figure, place plot hold on so multiple
    ## circles
    ## can be drawn
    # figure()
    clf
    hold on
    box on
    grid off

    mag_db = cellfun (@mag2db, mag, "uniformoutput", false);

    plot_args = horzcat (cellfun (@horzcat, pha, mag_db, sty, "uniformoutput", false){:});

    for k = 1:numsys
      numpoints = idivide (linspace (length ( mag_db{k} ) / 4, 3 * ...
      length ( mag_db{k} ) / 4, 10), int32 ( 1 ), "fix");
      for j = numpoints
        wstring = ["\\omega=", num2str(w{k}(j),"%.2f")];
        text (pha{k}(j), mag_db{k}(j), wstring, "rotation", 0);
      endfor
    endfor

    h11 = plot (plot_args{:});

    ## Nicer range and ticks for the phase
    x_lim = xlim ();
    y_lim = ylim ();
    if (rem (x_lim(1),30) != 0)
      x_lim(1) = x_lim(1) - (30 + rem (x_lim(1),30));
    endif
    if (x_lim(2) - x_lim(1) < 420)
      dx = 30;
    else
      dx = 60;
    endif
    xticks (x_lim(1):dx:x_lim(2));

    ## Now construct the grid with isobars for gain and phase
    ## of the closed loop

    col_pha = [0.6,0.6,0.6];
    col_mag = [0.6,0.6,0.6];
    pat_pha = ":";
    pat_mag = "-.";
    htxt = text (0,0,"-180°");
    ext = get (htxt, "extent");
    width_pha_label = ext(3);
    delete (htxt);

    ## Gain grid in dB

    ## Get min and max dB values for the grid
    db_min = y_lim(1) + sign (y_lim(1)) * rem (y_lim(1),20);
    dbgain = db_min:20:y_lim(2);
    db_default = [-10, -6, -3, -2, -1];
    dbgain = sort ([dbgain, db_default, -db_default, 0.5, 0.25, 0.1, 0.05]);

    ## Convert to true Gains
    M = 10 .^ ( dbgain ./ 20 );
    ## create an angle vector to draw the circles
    numsegments = 1024;
    omega = linspace ( -pi, pi ,numsegments ); ## numsegments segment circle
    ## important to go from -pi to pi to prevent wrap around lines on chart.
    for i = 1 : length ( dbgain );
      ## radius of the M-circles
      radius = M(i) / ( M(i) ^2 - 1 );
      ## Generate Plot Data for M-Circle with Real Portion on Y-axis
      Xm = -M(i) ^ 2 / ( M(i) ^ 2 - 1 ) + radius * cos ( omega );
      Ym = radius * sin( omega );
      ## Create polar plot (Phase Angle vs dB Gain)
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

      if dbgain(i) < 0
        tp1 = x_lim(2) - (x_lim(2)-x_lim(1))/100;
        pm = "-";
        ha = "right";
        v_off = 2;
      else
        pm = "+";
        tp1 = -180;
        ha = "center";
        v_off = 1;
      endif

      tp2 = db_plot_gain(1) + v_off;
      text (tp1, tp2, [pm,num2str(dbgain(i)),"dB"],...
            "color", col_mag, "clipping", "on",
            "horizontalalignment", ha);

      ## Plot the dB grid over the complete phase range
      x_min = x_lim(1);
      while (x_min < 0)
        plot ( phase, db_plot_gain, pat_mag, "color", col_mag, "linewidth", 0.4);
        x_min = x_min + 360;
        phase = phase - 360;
      endwhile

    endfor

    ## Phase grid in deg

    clphase = [ -178, -175, -170, -160, -150, -120, -90, -60, -30, -15, -10, -5, -2 ];
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

      ## Prepare labels for the phase
      if clphase(i) < 0
        pm = "-";
      else
        pm = "+";
      endif
      clphase_labels = clphase;

      ## Extend to lower nound of gein range
      [~,idx] = min (phase);
      phase = [phase(1:idx-2), clphase(i), clphase(i)-180, phase(idx+1:end)];
      db_plot_gain = [db_plot_gain(1:idx-2), y_lim(1)-1.5, y_lim(1)-1, db_plot_gain(idx+1:end)];

      ## Plot the phase grid over the complete phase range
      x_min = x_lim(1);
      while (x_min < 0)
        plot(phase, db_plot_gain, pat_pha, "color", col_pha, "linewidth", 0.4);
        x_min = x_min + 360;
        phase = phase - 360;
        for v = 1:2
          # Lables for all possible (periodic) occurrences
          text (clphase_labels(i)-1, y_lim(1) + 2*width_pha_label,...
                [pm,num2str(clphase_labels(i)),"°"],...
                "color", col_pha, "rotation", 90, "clipping", "on",...
                "horizontalalignment", "right");
          clphase_labels = clphase_labels - 180;
        endfor
      endwhile

    endfor

    ## Set correct ranges, title, labels and legend
    xlim (x_lim);
    ylim (y_lim);

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
