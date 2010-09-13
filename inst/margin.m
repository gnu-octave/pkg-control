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
## @deftypefn{Function File} {[@var{gamma}, @var{phi}, @var{w_gamma}, @var{w_phi}] =} margin (@var{sys})
## @deftypefnx{Function File} {[@var{gamma}, @var{phi}, @var{w_gamma}, @var{w_phi}] =} margin (@var{sys}, @var{tol})
## Gain and phase margin of a system. If no output arguments are given, both gain and phase margin
## are plotted on a bode diagram. Otherwise, the margins and their corresponding frequencies are
## computed and returned.
##
## @strong{Inputs}
## @table @var
## @item sys
## LTI model. Must be a single-input and single-output (SISO) system.
## @item tol
## Imaginary parts below @var{tol} are assumed to be zero.
## If not specified, default value @code{sqrt (eps)} is taken.
## @end table
##
## @strong{Outputs}
## @table @var
## @item gamma
## Gain margin (as gain, not dBs).
## @item phi
## Phase margin (in degrees).
## @item w_gamma
## Frequency for the gain margin (in rad/s).
## @item w_phi
## Frequency for the phase margin (in rad/s).
## @end table
##
## @seealso{roots}
##
## @example
## @group
## CONTINUOUS SYSTEMS
## Gain Margin
##         _               _
## L(jw) = L(jw)      BTW: L(jw) = L(-jw) = conj (L(jw))
##
## num(jw)   num(-jw)
## ------- = --------
## den(jw)   den(-jw)
##
## num(jw) den(-jw) = num(-jw) den(jw)
##
## imag (num(jw) den(-jw)) = 0
## imag (num(-jw) den(jw)) = 0
##
## Phase Margin
##           |num(jw)|
## |L(jw)| = |-------| = 1
##           |den(jw)|
##   _     2      2
## z z = Re z + Im z
##
## num(jw)   num(-jw)
## ------- * -------- = 1
## den(jw)   den(-jw)
##
## num(jw) num(-jw) - den(jw) den(-jw) = 0
##
## real (num(jw) num(-jw) - den(jw) den(-jw)) = 0
##
##
## DISCRETE SYSTEMS
## Gain Margin
##                              jwT         log z
## L(z) = L(1/z)      BTW: z = e    --> w = -----
##                                           j T
## num(z)   num(1/z)
## ------ = --------
## den(z)   den(1/z)
##
## num(z) den(1/z) - num(1/z) den(z) = 0
##
## Phase Margin
##          |num(z)|
## |L(z)| = |------| = 1
##          |den(z)|
##
## L(z) L(1/z) = 1
##
## num(z)   num(1/z)
## ------ * -------- = 1
## den(z)   den(1/z)
##
## num(z) num(1/z) - den(z) den(1/z) = 0
##
## PS: How to get L(1/z)
##           4       3       2
## p(z) = a z  +  b z  +  c z  +  d z  +  e
##
##             -4      -3      -2      -1
## p(1/z) = a z  +  b z  +  c z  +  d z  +  e
##
##           -4                    2       3       4
##        = z   ( a  +  b z  +  c z  +  d z  +  e z  )
##
##               4       3       2                     4
##        = ( e z  +  d z  +  c z  +  b z  +  a ) / ( z  )
## @end group
## @end example
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: July 2009
## Version: 0.7.1

function [gamma_r, phi_r, w_gamma_r, w_phi_r] = margin (sys, tol = sqrt (eps))

  ## check whether arguments are OK
  if (nargin < 1 || nargin > 2)
    print_usage ();
  endif

  if (! isa (sys, "lti"))
    error ("margin: argument sys must be a LTI system");
  endif

  if (! issiso (sys))
    error ("margin: argument sys must be a SISO system");
  endif

  ## get transfer function
  [num, den, Ts] = tfdata (sys);
  num = num{1, 1};
  den = den{1, 1};


  if (Ts == 0)  # CONTINUOUS SYSTEM

    ## create polynomials s -> jw
    l_num = length (num);
    l_den = length (den);
    num_jw = zeros (1, l_num);
    den_jw = zeros (1, l_den);

    for k = 1 : l_num
      num_jw(k) = num(k) * i^(l_num - k);
    endfor

    for k = 1 : l_den
      den_jw(k) = den(k) * i^(l_den - k);
    endfor

    ## GAIN MARGIN
    ## create gm polynomial
    gm_poly = imag (conv (num_jw, conj (den_jw)));

    ## find frequencies w
    w = roots (gm_poly);

    ## filter results
    [gamma, w_gamma] = gm_filter (w, num, den, Ts, tol);

    ## PHASE MARGIN
    ## create pm polynomials
    poly_1 = conv (num_jw, conj (num_jw));
    poly_2 = conv (den_jw, conj (den_jw));

    ## make polynomials equally long for subtraction
    [poly_1, poly_2] = poly_equalizer (poly_1, poly_2);

    ## subtract polynomials
    pm_poly = real (poly_1 - poly_2);

    ## find frequencies w
    w = roots (pm_poly);

    ## filter results
    [phi, w_phi] = pm_filter (w, num, den, Ts, tol);


  else  # DISCRETE SYSTEM

    ## create polynomials z -> 1/z
    l_num = length (num);
    l_den = length (den);

    num_rev = fliplr (num);
    den_rev = fliplr (den);

    num_div = zeros (1, l_num);
    den_div = zeros (1, l_den);
    num_div(1) = 1;
    den_div(1) = 1;

    num_inv = conv (num_rev, den_div);
    den_inv = conv (den_rev, num_div);

    ## GAIN MARGIN
    ## create gm polynomial
    poly_1 = conv (num, den_inv);
    poly_2 = conv (num_inv, den);

    ## make polynomials equally long for subtraction
    [poly_1, poly_2] = poly_equalizer (poly_1, poly_2);

    ## subtract polynomials
    gm_poly = poly_1 - poly_2;

    ## find frequencies z
    z = roots (gm_poly);

    ## filter results
    idx = find (abs (abs (z) - 1) < tol);  # find z with magnitude 1

    if (length (idx) > 0)  # if z with magnitude 1 exist
      z_gm = z(idx);
      w = log (z_gm) / (i*Ts);  # get frequencies w from z
      [gamma, w_gamma] = gm_filter (w, num, den, Ts, tol);
    else  # there are no z with magnitude 1
      gamma = Inf;
      w_gamma = NaN;
    endif

    ## PHASE MARGIN
    ## create pm polynomials
    poly_1 = conv (num, num_inv);
    poly_2 = conv (den, den_inv);

    ## make polynomials equally long for subtraction
    [poly_1, poly_2] = poly_equalizer (poly_1, poly_2);

    ## subtract polynomials
    pm_poly = poly_1 - poly_2;

    ## find frequencies z
    z = roots (pm_poly);

    ## filter results
    idx = find (abs (abs (z) - 1) < tol);  # find z with magnitude 1

    if (length (idx) > 0)  # if z with magnitude 1 exist
      z_gm = z(idx);
      w = log (z_gm) / (i*Ts);  # get frequencies w from z
      [phi, w_phi] = pm_filter (w, num, den, Ts, tol);
    else  # there are no z with magnitude 1
      phi = 180;
      w_phi = NaN;
    endif

  endif


  if (nargout == 0)  # show bode diagram

    [H, w] = __getfreqresp__ (sys, [], false, 0, "std");

    H = reshape (H, [], 1);
    mag_db = 20 * log10 (abs (H));
    pha = unwrap (arg (H)) * 180 / pi;
    gamma_db = 20 * log10 (gamma);

    wv = [min(w), max(w)];
    ax_vec_mag = __axis2dlim__ ([reshape(w, [], 1), reshape(mag_db, [], 1)]);
    ax_vec_mag(1:2) = wv;
    ax_vec_pha = __axis2dlim__ ([reshape(w, [], 1), reshape(pha, [], 1)]);
    ax_vec_pha(1:2) = wv;

    wgm = [w_gamma, w_gamma];
    mgmh = [-gamma_db, ax_vec_mag(3)];
    mgm = [0, -gamma_db];
    pgm = [ax_vec_pha(4), -180];

    wpm = [w_phi, w_phi];
    mpm = [0, ax_vec_mag(3)];
    ppmh = [ax_vec_pha(4), phi - 180];
    ppm = [phi - 180, -180];

    title_str = sprintf ("GM = %g dB (at %g rad/s),   PM = %g deg (at %g rad/s)",
                         gamma_db, w_gamma, phi, w_phi);
    if (Ts == 0)
      xl_str = "Frequency [rad/s]";
    else
      xl_str = sprintf ("Frequency [rad/s]     w_N = %g", pi/Ts);
    endif

    subplot (2, 1, 1)
    semilogx (w, mag_db, "b", wv, [0, 0], ":k", wgm, mgmh, ":k", wgm, mgm, "r", wpm, mpm, ":k")
    axis (ax_vec_mag)
    grid ("on")
    title (title_str)
    ylabel ("Magnitude [dB]")

    subplot (2, 1, 2)
    semilogx (w, pha, "b", wv, [-180, -180], ":k", wgm, pgm, ":k", wpm, ppmh, ":k", wpm, ppm, "r")
    axis (ax_vec_pha)
    grid ("on")
    xlabel (xl_str)
    ylabel ("Phase [deg]")

  else  # return values
    gamma_r = gamma;
    phi_r = phi;
    w_gamma_r = w_gamma;
    w_phi_r = w_phi;
  endif

endfunction


function [poly_eq_1, poly_eq_2] = poly_equalizer (poly_1, poly_2)

  l_p1 = length (poly_1);
  l_p2 = length (poly_2);
  l_max = max (l_p1, l_p2);

  lead_zer_1 = zeros (1, l_max - l_p1);
  lead_zer_2 = zeros (1, l_max - l_p2);

  poly_eq_1 = horzcat (lead_zer_1, poly_1);
  poly_eq_2 = horzcat (lead_zer_2, poly_2);

endfunction


function [gamma, w_gamma] = gm_filter (w, num, den, Ts, tol)

  idx = find ((abs (imag (w)) < tol) & (real (w) > 0));  # find frequencies in R+
  l_idx = length (idx);

  if (l_idx > 0)  # if frequencies in R+ exist
    w_gm = real (w(idx));
    f_resp = zeros (1, l_idx);
    gm = zeros (1, l_idx);

    if (Ts == 0)
      for k = 1 : l_idx
        f_resp(k) = polyval (num, i*w_gm(k)) / polyval (den, i*w_gm(k));
        gm(k) = inv (abs (f_resp(k)));
      endfor
    else
      for k = 1 : l_idx
        f_resp(k) = polyval (num, exp (i*w_gm(k)*Ts)) / polyval (den, exp (i*w_gm(k)*Ts));
        gm(k) = inv (abs (f_resp(k)));
      endfor
    endif

    ## find crossings between 0 and -1
    idx = find ((real (f_resp) < 0) & (real (f_resp) >= -1));

    if (length (idx) > 0)  # if crossings between 0 and -1 exist
      gm = gm(idx);
      w_gm = w_gm(idx);
      [gamma, idx] = min (gm);
      w_gamma = w_gm(idx);
    else  # there are no crossings between 0 and -1
      idx = find (real (f_resp) < -1);  # find crossings between -1 and -Inf

      if (length (idx) > 0)  # if crossings between -1 and -Inf exist
        gm = gm(idx);
        w_gm = w_gm(idx);
        [gamma, idx] = max (gm);
        w_gamma = w_gm(idx);
      else
        gamma = Inf;
        w_gamma = NaN;
      endif
    endif
  else  # there are no frequencies in R+
    gamma = Inf;
    w_gamma = NaN;
  endif

endfunction


function [phi, w_phi] = pm_filter (w, num, den, Ts, tol)

  idx = find ((abs (imag (w)) < tol) & (real (w) > 0));  # find frequencies in R+
  l_idx = length (idx);

  if (l_idx > 0)  # if frequencies in R+ exist
    w_pm = real (w(idx));
    pm = zeros (1, l_idx);

    if (Ts == 0)
      for k = 1 : l_idx
        f_resp = polyval (num, i*w_pm(k)) / polyval (den, i*w_pm(k));
        pm(k) = 180  +  arg (f_resp) / pi * 180;
      endfor
    else
      for k = 1 : l_idx
        f_resp = polyval (num, exp (i*w_pm(k)*Ts)) / polyval (den, exp (i*w_pm(k)*Ts));
        pm(k) = 180  +  arg (f_resp) / pi * 180;
      endfor
    endif

    [phi, idx] = min (pm);
    w_phi = w_pm(idx);
  else  # there are no frequencies in R+
    phi = 180;
    w_phi = NaN;
  endif

endfunction


%!shared margin_c, margin_c_exp, margin_d, margin_d_exp
%! sysc = tf ([24], [1, 6, 11, 6]);
%! [gamma_c, phi_c, w_gamma_c, w_phi_c] = margin (sysc);
%! sysd = tf (c2d (ss (sysc), 0.3));  # c2d for tf not implemented yet
%! [gamma_d, phi_d, w_gamma_d, w_phi_d] = margin (sysd);
%!
%! margin_c = [gamma_c, phi_c, w_gamma_c, w_phi_c];
%! margin_d = [gamma_d, phi_d, w_gamma_d, w_phi_d];
%!
%! ## results from this implementation and the "dark side" diverge
%! ## from the third digit after the decimal point on
%!
%! gamma_c_exp = 2.50;
%! phi_c_exp = 35.43;
%! w_gamma_c_exp = 3.32;
%! w_phi_c_exp = 2.06;
%!
%! gamma_d_exp = 1.41;
%! phi_d_exp = 18.60;
%! w_gamma_d_exp = 2.48;
%! w_phi_d_exp = 2.04;
%!
%! margin_c_exp = [gamma_c_exp, phi_c_exp, w_gamma_c_exp, w_phi_c_exp];
%! margin_d_exp = [gamma_d_exp, phi_d_exp, w_gamma_d_exp, w_phi_d_exp];
%!
%!assert (margin_c, margin_c_exp, 1e-2);
%!assert (margin_d, margin_d_exp, 1e-2);
