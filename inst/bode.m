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
## the response is printed on the screen. In the latter case, the chidlren
## of the figure handle are the handles for (1) the phase plot, (2) the
## legend, and (3) the magnitude plot. The magnitude plot is the active
## handle after the bode plot is finished.
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

  names = arrayfun (@inputname, 1:nargin, 'uniformoutput', false);

  [H, w, sty, sys_idx, H_auto, w_auto, leg] = __frequency_response__ ("bode", varargin, nargout, names);

  H = cellfun (@reshape, H, {[]}, {1}, "uniformoutput", false);
  H_auto = cellfun (@reshape, H_auto, {[]}, {1}, "uniformoutput", false);
  mag = cellfun (@abs, H, "uniformoutput", false);
  pha = cellfun (@(H) unwrap (arg (H)) * 180 / pi, H, "uniformoutput", false);
  pha_auto = cellfun (@(H_auto) unwrap (arg (H_auto)) * 180 / pi, H_auto, "uniformoutput", false);
  numsys = length (sys_idx);

  ## check for poles and zeroes at the origin (or at one for discrete-time
  ## system) for each of the numsys systems and compare this asymptotic
  ## value (omega -> 0) to initial phase for the full auto w-range, which
  ## is supposed to start at sufficiently small values fo w
  initial_phase_auto = cellfun(@(x) x(1), pha_auto);
  initial_phase = cellfun(@(x) x(1), pha);

  all_poles = cell(1,numsys);
  all_zeros = cell(1,numsys);

  for h = 1:numsys

    sys1 = varargin{sys_idx(h)};
    [num,den] = tfdata (sys1,'v');
    [all_zeros(h), all_poles(h), k] = zpkdata (sys1);

    pole_integrator = 0;
    if (isdt (sys1))
      pole_integrator = 1;
    endif
    tol = sqrt (eps);

    n_poles_at_origin = sum (abs (all_poles{h} - pole_integrator) < tol);
    n_zeros_at_origin = sum (abs (all_zeros{h} - pole_integrator) < tol);
    asymptotic_low_freq_phase = (n_zeros_at_origin - n_poles_at_origin)*90;
    if (k < 0)
      asymptotic_low_freq_phase = asymptotic_low_freq_phase - 180;
    endif

    pha_auto(h)={cell2mat(pha_auto(h)) + round((asymptotic_low_freq_phase - initial_phase_auto(h))/90)*90};
    pha(h)={cell2mat(pha(h)) + round((asymptotic_low_freq_phase - initial_phase(h))/90)*90};

  endfor

  ## check if custom w-range is at higher frequencies and provide and if yes,
  ## add missing "history" for the phase
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
        pha(h)={cell2mat(pha(h)) + round((pha_cmp - initial_phase(h))/90)*90};
      endif
    else
      ## w range is identical to auto range, just use the phase for the
      ## auto range
      pha(h) = pha_auto(h);
    endif
  endfor

  ## now fix the phase in case of a conjugate complex pair of poles on the
  ## imag. axis

  tol = sqrt (eps);
  max_idx = length (pha{1});

  w0_poles = cellfun (@(p) imag (p(abs(real(p))<tol & imag(p)>0)), all_poles,...
                      'uniformoutput', false);
  w0_zeros = cellfun (@(p) imag (p(abs(real(p))<tol & imag(p)>0)), all_zeros,...
                      'uniformoutput', false);

  dpha = cellfun (@diff, pha, 'uniformoutput', false);

  w0_inc = cellfun (@(dp) find (dp > 90), dpha, 'uniformoutput', false);
  w0_dec = cellfun (@(dp) find (dp < -90), dpha, 'uniformoutput', false);


  for h = 1:numsys

    for j = 1:length(w0_inc{h})
      if (w0_inc{h}(j)-2 > 0 && w0_inc{h}(j)+2 <= max_idx)
        lb = w{h}(w0_inc{h}(j)-2) < w0_poles{h};
        ub = w{h}(w0_inc{h}(j)+2) > w0_poles{h};
        if sum (lb & ub) == 1
          ## phase is increasing although poles
          pha{h}(w0_inc{h}(j)+1:end) -= 360;
        endif
      endif
    endfor

    for j = 1:length(w0_dec{h})
      if (w0_dec{h}(j)-2 > 0 && w0_dec{h}(j)+2 <= max_idx)
        lb = w{h}(w0_dec{h}(j)-2) < w0_zeros{h};
        ub = w{h}(w0_dec{h}(j)+2) > w0_zeros{h};
        if sum (lb & ub) == 1
          ## phase is decreasing although zeros
          pha{h}(w0_dec{h}(j)+1:end) += 360;
        endif
      endif
    endfor

  endfor


  ## do the plotting

  if (! nargout)

    ## plot
    mag_db = cellfun (@mag2db, mag, "uniformoutput", false);

    mag_args = horzcat (cellfun (@horzcat, w, mag_db, sty, "uniformoutput", false){:});
    pha_args = horzcat (cellfun (@horzcat, w, pha, sty, "uniformoutput", false){:});

    hax1 = subplot (2, 1, 1);
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
    legend ("off")

    ## Make top axes current for ML compatibility
    set (gcf, "currentaxes", hax1);

  else

    ## no plotting, assign values to the output parameters
    mag_r = mag{1};
    pha_r = pha{1};
    w_r = w{1};

  endif

endfunction

%!demo
%! G_1 = tf ([10 1],[1 0.1 1 0]);
%! G_2 = tf ([1 0.1 1],[100 1 1]);
%! bode(G_1,G_2);

%!test
%! s = tf('s');
%! K = 1;
%! T = 2;
%! G = K/(1+T*s);
%! [mag phas w] = bode(G);
%! mag_dB = 20*log10(mag);
%! index = find(mag_dB < -3,1);
%! w_cutoff = w(index);
%! assert (1/T, w_cutoff, eps);

%!test
%! G = tf ([100 0 1],conv([1 0.00001 1],[1/100 0 1]));
%! [mag pha w] = bode(G);
%! pha_final = pha(end);
%! assert (-180, pha_final, 1e-4);

%!test
%! s = tf ('s');
%! G = s/(1+s)^4;
%! # every whole decade has an w in the list
%! w = logspace (-4,2,61);
%! [mag, pha, w] = bode (G, w);
%! wE = 1;
%! pha_w0 = atan2 (1e-4,4e-8) / pi * 180;
%! abs_w0 = 1e-4 / sqrt (1^2 + 1e-4^2)^4;
%! pha_wE = -90;
%! abs_wE = 1/4;
%! assert (pha(1),  pha_w0, 1e-8);
%! assert (mag(1),  abs_w0, 1e-8);
%! assert (pha(41), pha_wE, 1e-8);
%! assert (mag(41), abs_wE, 1e-8);

