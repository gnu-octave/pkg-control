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
##           |num(jw)|        _
## |L(jw)| = |-------| = 1 = |L(jw)|
##           |den(jw)|
##   _     2      2
## z z = Re z + Im z
##
## num(jw) num(-jw)
## ---------------- = 1
## den(jw) den(-jw)
##
## num(jw) num(-jw) - den(jw) den(-jw) = 0
##
## real (num(jw) num(-jw) - den(jw) den(-jw)) = 0
## @end group
## @end example
## @end deftypefn

## Version: 0.3

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
    idx_1 = find (abs (imag (w)) < tol);  # find real frequencies
    idx_2 = find (real (w) > 0);  # find positive frequencies
    idx = intersect (idx_1, idx_2);  # find frequencies in R+

    if (length (idx) > 0)  # if frequencies in R+ exist
      w_gm = real (w(idx));

      for k = 1 : length (w_gm)
        f_resp(k) = polyval (num_jw, w_gm(k)) / polyval (den_jw, w_gm(k));
        gm(k) = inv (abs (f_resp(k)));
      endfor

      ## find crossings between 0 and -1
      idx_3 = find (real (f_resp) < 0);
      idx_4 = find (real (f_resp) >= -1);
      idx_5 = intersect (idx_3, idx_4);

      if (length (idx_5) > 0)  # if crossings between 0 and -1 exist
        gm = gm(idx_5);
        w_gm = w_gm(idx_5);

        [gamma, idx_6] = min (gm);
        w_gamma = w_gm(idx_6);
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

    poly_eq_1 = zeros (1, l_max);
    poly_eq_2 = zeros (1, l_max);

    for k = 1 : l_max
      if (l_max - k < l_p1)
        poly_eq_1(k) = poly_1(k + l_p1 - l_max);
      endif

      if (l_max - k < l_p2)
        poly_eq_2(k) = poly_2(k + l_p2 - l_max);
      endif
    endfor

    ## subtract polynomials
    pm_poly = real (poly_eq_1 - poly_eq_2);

    ## find frequencies w
    w = roots (pm_poly);

    ## filter results
    idx_1 = find (abs (imag (w)) < tol);  # find real frequencies
    idx_2 = find (real (w) > 0);  # find positive frequencies
    idx = intersect (idx_1, idx_2);  # find frequencies in R+

    if (length (idx) > 0)  # if frequencies in R+ exist
      w_pm = real (w(idx));

      for k = 1 : length (w_pm)
        f_resp = polyval (num_jw, w_pm(k)) / polyval (den_jw, w_pm(k));
        pm(k) = 180  +  arg (f_resp) / pi * 180;
      endfor

      [phi, idx_3] = min (pm);
      w_phi = w_pm(idx_3);
    else  # there are no frequencies in R+
      phi = 180;
      w_phi = NaN;
    endif


  else  # DISCRETE SYSTEM

    ## create polynomials z -> 1/z
    l_num = length (num);
    l_den = length (den);
    num_inv = zeros (1, l_num);
    den_inv = zeros (1, l_den);

    for k = 1 : l_num
      num_inv(k) = num(l_num + 1 - k);
    endfor

    for k = 1 : l_den
      den_inv(k) = den(l_den + 1 - k);
    endfor

    z_num = zeros (1, l_num);
    z_den = zeros (1, l_den);
    z_num(1) = 1;
    z_den(1) = 1; 

    ## GAIN MARGIN
    ## create gm polynomial
    poly_1 = conv (conv (num, den_inv), z_num);
    poly_2 = conv (conv (num_inv, den), z_den);

    ## make polynomials equally long for subtraction
    l_p1 = length (poly_1);
    l_p2 = length (poly_2);
    l_max = max (l_p1, l_p2);

    poly_eq_1 = zeros (1, l_max);
    poly_eq_2 = zeros (1, l_max);

    for k = 1 : l_max
      if (l_max - k < l_p1)
        poly_eq_1(k) = poly_1(k + l_p1 - l_max);
      endif

      if (l_max - k < l_p2)
        poly_eq_2(k) = poly_2(k + l_p2 - l_max);
      endif
    endfor

    ## subtract polynomials
    gm_poly = poly_eq_1 - poly_eq_2;

    ## find frequencies z
    z = roots (gm_poly);

    ## filter results
    idx = find (abs (abs (z) - 1) < tol);  # find z with magnitude 1

    if (length (idx) > 0)  # if z with magnitude 1 exist
      z_gm = z(idx);
      w = log (z_gm) / (i*Ts);  # get frequencies w from z

      idx_1 = find (abs (imag (w)) < tol);  # find real frequencies
      idx_2 = find (real (w) > 0);  # find positive frequencies
      idx = intersect (idx_1, idx_2);  # find frequencies in R+

      if (length (idx) > 0)  # if frequencies in R+ exist
        w_gm = real (w(idx));

        for k = 1 : length (w_gm)
          f_resp(k) = polyval (num, exp (i*w_gm(k)*Ts)) / polyval (den, exp (i*w_gm(k)*Ts));
          gm(k) = inv (abs (f_resp(k)));
        endfor

        ## find crossings between 0 and -1
        idx_3 = find (real (f_resp) < 0);
        idx_4 = find (real (f_resp) >= -1);
        idx_5 = intersect (idx_3, idx_4);

        if (length (idx_5) > 0)  # if crossings between 0 and -1 exist
          gm = gm(idx_5);
          w_gm = w_gm(idx_5);

          [gamma, idx_6] = min (gm);
          w_gamma = w_gm(idx_6);
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
    poly_1 = conv (conv (num, num_inv), z_den);
    poly_2 = conv (conv (den, den_inv), z_num);

    ## make polynomials equally long for subtraction
    l_p1 = length (poly_1);
    l_p2 = length (poly_2);
    l_max = max (l_p1, l_p2);

    poly_eq_1 = zeros (1, l_max);
    poly_eq_2 = zeros (1, l_max);

    for k = 1 : l_max
      if (l_max - k < l_p1)
        poly_eq_1(k) = poly_1(k + l_p1 - l_max);
      endif

      if (l_max - k < l_p2)
        poly_eq_2(k) = poly_2(k + l_p2 - l_max);
      endif
    endfor

    ## subtract polynomials
    pm_poly = poly_eq_1 - poly_eq_2;

    ## find frequencies z
    z = roots (pm_poly);

    ## filter results
    idx = find (abs (abs (z) - 1) < tol);  # find z with magnitude 1

    if (length (idx) > 0)  # if z with magnitude 1 exist
      z_gm = z(idx);
      w = log (z_gm) / (i*Ts);  # get frequencies w from z

      idx_1 = find (abs (imag (w)) < tol);  # find real frequencies
      idx_2 = find (real (w) > 0);  # find positive frequencies
      idx = intersect (idx_1, idx_2);  # find frequencies in R+

      if (length (idx) > 0)  # if frequencies in R+ exist
        w_pm = real (w(idx));

        for k = 1 : length (w_pm)
          f_resp = polyval (num, exp (i*w_pm(k)*Ts)) / polyval (den, exp (i*w_pm(k)*Ts));
          pm(k) = 180  +  arg (f_resp) / pi * 180;
        endfor

        [phi, idx_3] = min (pm);
        w_phi = w_pm(idx_3);
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


%!shared gamma, phi, w_gamma, w_phi, gamma_exp, phi_exp, w_gamma_exp, w_phi_exp
%! sys = tf ([24], [1, 6, 11, 6]);
%! [gamma, phi, w_gamma, w_phi] = margin (sys);
%!
%! gamma_exp = 2.50000000000000;
%! phi_exp = 35.4254199887541;
%! w_gamma_exp = 3.31662479035540;
%! w_phi_exp = 2.06385145800825;
%!
%!assert (gamma, gamma_exp, 10*eps);
%!assert (phi, phi_exp, 256*eps); # FIXME: why this high tol?
%!assert (w_gamma, w_gamma_exp, 2*eps);
%!assert (w_phi, w_phi_exp, 25*eps);
