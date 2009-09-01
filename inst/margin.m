## Copyright (C) 2009 Lukas Reichlin <lukas.reichlin@swissonline.ch>
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn{Function File} {[@var{gamma}, @var{phi}, @var{w_gamma}, @var{w_phi}] =} margin (@var{sys})
## @deftypefnx{Function File} {[@var{gamma}, @var{phi}, @var{w_gamma}, @var{w_phi}] =} margin (@var{sys}, @var{tol})
## Gain and phase margin of a system.
##
## @strong{Inputs}
## @table @var
## @item sys
## System data structure of a SISO system.
## @item tol
## Imaginary parts below tol are assumed to be zero. If not specified, default value 1e-7 is taken.
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

## Version: 0.4.4

function [gamma, phi, w_gamma, w_phi] = margin (sys, tol)

  ## check whether arguments are OK
  if (nargin < 1 || nargin > 2)
    print_usage ();
  endif

  if (! isstruct (sys))
    error ("margin: argument sys must be a system data structure");
  endif

  if (! is_siso (sys))
    error ("margin: argument sys must be a SISO system");
  endif

  if (is_digital (sys, 2) == -1)
    error ("margin: system must be either purely continuous or purely discrete");
  endif

  if (nargin == 1)
    tol = 1e-7;
  endif

  ## get transfer function
  [num, den, Ts] = sys2tf (sys);


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
    idx = find ((abs (imag (w)) < tol) & (real (w) > 0));  # find frequencies in R+
    l_idx = length (idx);

    if (l_idx > 0)  # if frequencies in R+ exist
      w_gm = real (w(idx));
      f_resp = zeros (1, l_idx);
      gm = zeros (1, l_idx);

      for k = 1 : l_idx
        f_resp(k) = polyval (num_jw, w_gm(k)) / polyval (den_jw, w_gm(k));
        gm(k) = inv (abs (f_resp(k)));
      endfor

      ## find crossings between 0 and -1
      idx = find ((real (f_resp) < 0) & (real (f_resp) >= -1));

      if (length (idx) > 0)  # if crossings between 0 and -1 exist
        gm = gm(idx);
        w_gm = w_gm(idx);

        [gamma, idx] = min (gm);
        w_gamma = w_gm(idx);
      else  # there are no crossings between 0 and -1
        gamma = Inf;
        w_gamma = NaN;
      endif
    else  # there are no frequencies in R+
      gamma = Inf;
      w_gamma = NaN;
    endif

    ## PHASE MARGIN
    ## create pm polynomials
    poly_1 = conv (num_jw, conj (num_jw));
    poly_2 = conv (den_jw, conj (den_jw));

    ## make polynomials equally long for subtraction
    l_p1 = length (poly_1);
    l_p2 = length (poly_2);
    l_max = max (l_p1, l_p2);

    lead_zer_1 = zeros (1, l_max - l_p1);
    lead_zer_2 = zeros (1, l_max - l_p2);

    poly_eq_1 = horzcat (lead_zer_1, poly_1);
    poly_eq_2 = horzcat (lead_zer_2, poly_2);

    ## subtract polynomials
    pm_poly = real (poly_eq_1 - poly_eq_2);

    ## find frequencies w
    w = roots (pm_poly);

    ## filter results
    idx = find ((abs (imag (w)) < tol) & (real (w) > 0));  # find frequencies in R+
    l_idx = length (idx);

    if (l_idx > 0)  # if frequencies in R+ exist
      w_pm = real (w(idx));
      pm = zeros (1, l_idx);

      for k = 1 : l_idx
        f_resp = polyval (num_jw, w_pm(k)) / polyval (den_jw, w_pm(k));
        pm(k) = 180  +  arg (f_resp) / pi * 180;
      endfor

      [phi, idx] = min (pm);
      w_phi = w_pm(idx);
    else  # there are no frequencies in R+
      phi = 180;
      w_phi = NaN;
    endif


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
    l_p1 = length (poly_1);
    l_p2 = length (poly_2);
    l_max = max (l_p1, l_p2);

    lead_zer_1 = zeros (1, l_max - l_p1);
    lead_zer_2 = zeros (1, l_max - l_p2);

    poly_eq_1 = horzcat (lead_zer_1, poly_1);
    poly_eq_2 = horzcat (lead_zer_2, poly_2);

    ## subtract polynomials
    gm_poly = poly_eq_1 - poly_eq_2;

    ## find frequencies z
    z = roots (gm_poly);

    ## filter results
    idx = find (abs (abs (z) - 1) < tol);  # find z with magnitude 1

    if (length (idx) > 0)  # if z with magnitude 1 exist
      z_gm = z(idx);
      w = log (z_gm) / (i*Ts);  # get frequencies w from z

      idx = find ((abs (imag (w)) < tol) & (real (w) > 0));  # find frequencies in R+
      l_idx = length (idx);

      if (l_idx > 0)  # if frequencies in R+ exist
        w_gm = real (w(idx));
        f_resp = zeros (1, l_idx);
        gm = zeros (1, l_idx);

        for k = 1 : l_idx
          f_resp(k) = polyval (num, exp (i*w_gm(k)*Ts)) / polyval (den, exp (i*w_gm(k)*Ts));
          gm(k) = inv (abs (f_resp(k)));
        endfor

        ## find crossings between 0 and -1
        idx = find ((real (f_resp) < 0) & (real (f_resp) >= -1));

        if (length (idx) > 0)  # if crossings between 0 and -1 exist
          gm = gm(idx);
          w_gm = w_gm(idx);

          [gamma, idx] = min (gm);
          w_gamma = w_gm(idx);
        else  # there are no crossings between 0 and -1
          gamma = Inf;
          w_gamma = NaN;
        endif
      else  # there are no frequencies in R+
        gamma = Inf;
        w_gamma = NaN;
      endif
    else  # there are no z with magnitude 1
      gamma = Inf;
      w_gamma = NaN;
    endif

    ## PHASE MARGIN
    ## create pm polynomials
    poly_1 = conv (num, num_inv);
    poly_2 = conv (den, den_inv);

    ## make polynomials equally long for subtraction
    l_p1 = length (poly_1);
    l_p2 = length (poly_2);
    l_max = max (l_p1, l_p2);

    lead_zer_1 = zeros (1, l_max - l_p1);
    lead_zer_2 = zeros (1, l_max - l_p2);

    poly_eq_1 = horzcat (lead_zer_1, poly_1);
    poly_eq_2 = horzcat (lead_zer_2, poly_2);

    ## subtract polynomials
    pm_poly = poly_eq_1 - poly_eq_2;

    ## find frequencies z
    z = roots (pm_poly);

    ## filter results
    idx = find (abs (abs (z) - 1) < tol);  # find z with magnitude 1

    if (length (idx) > 0)  # if z with magnitude 1 exist
      z_gm = z(idx);
      w = log (z_gm) / (i*Ts);  # get frequencies w from z

      idx = find ((abs (imag (w)) < tol) & (real (w) > 0));  # find frequencies in R+
      l_idx = length (idx);

      if (l_idx > 0)  # if frequencies in R+ exist
        w_pm = real (w(idx));
        pm = zeros (1, l_idx);

        for k = 1 : l_idx
          f_resp = polyval (num, exp (i*w_pm(k)*Ts)) / polyval (den, exp (i*w_pm(k)*Ts));
          pm(k) = 180  +  arg (f_resp) / pi * 180;
        endfor

        [phi, idx] = min (pm);
        w_phi = w_pm(idx);
      else  # there are no frequencies in R+
        phi = 180;
        w_phi = NaN;
      endif
    else  # there are no z with magnitude 1
      phi = 180;
      w_phi = NaN;
    endif

  endif

endfunction


%!shared margin_c, margin_c_exp, margin_d, margin_d_exp
%! sysc = tf ([24], [1, 6, 11, 6]);
%! [gamma_c, phi_c, w_gamma_c, w_phi_c] = margin (sysc);
%! sysd = c2d (sysc, 0.3);
%! [gamma_d, phi_d, w_gamma_d, w_phi_d] = margin (sysd);
%!
%! gamma_c_exp = 2.50000000000000;
%! phi_c_exp = 35.4254199887541;
%! w_gamma_c_exp = 3.31662479035540;
%! w_phi_c_exp = 2.06385145800825;
%!
%! gamma_d_exp = 1.41205400689315;
%! phi_d_exp = 18.6021379182828;
%! w_gamma_d_exp = 2.47750745165538;
%! w_phi_d_exp = 2.04369751003386;
%!
%! margin_c = [gamma_c, phi_c, w_gamma_c, w_phi_c];
%! margin_c_exp = [gamma_c_exp, phi_c_exp, w_gamma_c_exp, w_phi_c_exp];
%! margin_d = [gamma_d, phi_d, w_gamma_d, w_phi_d];
%! margin_d_exp = [gamma_d_exp, phi_d_exp, w_gamma_d_exp, w_phi_d_exp];
%!
%!assert (margin_c, margin_c_exp, 256*eps); # FIXME: why such a high tol?
%!assert (margin_d, margin_d_exp, 2976*eps); # FIXME: why such a high tol?
