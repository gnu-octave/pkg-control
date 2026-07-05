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
## Frequency response of TF models.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.5

function H = __freqresp__ (sys, w, cellflag = false)

  [num, den, tsam] = tfdata (sys, "vector");

  ## Only use tf representation for low order systems without
  ## extremely spreaded coefficients

  use_ss = true;

  if (issiso (sys)) && (length (den) < 6)

    spread_num = max (log10 (abs (num(abs(num)>0)))) ...
                 - min (log10 (abs (num(abs(num)>0))));
    spread_den = max (log10 (abs (den(abs(den)>0)))) ...
                 - min (log10 (abs (den(abs(den)>0))));

    if (spread_num < 8) && (spread_den < 8)
      use_ss = false;
    endif

  endif

  if (use_ss)

    if (any (w == 0))
      dc_idx = find (w(:) == 0);
      H_dc = __freqresp_tf__ (sys, num, den, tsam, 0, false);

      if (numel (dc_idx) == numel (w))
        if (cellflag)
          H = repmat ({H_dc}, numel (w), 1);
        else
          H = repmat (H_dc, [1, 1, numel(w)]);
        endif
      else
        nz_idx = find (w(:) != 0);
        sys_ss = ss (sys);
        H_nz = __freqresp__ (sys_ss, w(nz_idx), cellflag);

        if (cellflag)
          H = cell (size (w));
          H(dc_idx) = repmat ({H_dc}, numel (dc_idx), 1);
          H(nz_idx) = H_nz(:);
        else
          [p, m] = size (sys);
          H = zeros (p, m, numel (w));
          H(:, :, dc_idx) = repmat (H_dc, [1, 1, numel(dc_idx)]);
          H(:, :, nz_idx) = H_nz;
        endif
      endif

      return;
    endif

    sys_ss = ss (sys);
    H = __freqresp__ (sys_ss, w, cellflag);

    return;

  endif


  H = __freqresp_tf__ (sys, num, den, tsam, w, cellflag);

endfunction


function H = __freqresp_tf__ (sys, num, den, tsam, w, cellflag)

  ## Use tf represantation

  if (isct (sys))  # continuous system
    s = i * w;
  else             # discrete system
    s = exp (i * w * abs (tsam));
  endif

  s = reshape (s, 1, 1, []);

  if (issiso (sys))
    H = __fracval__ (num, den, s);
  else
    H = cellfun (@(x, y) __fracval__ (x, y, s), num, den, "uniformoutput", false);
    H = cell2mat (H);
  endif

  if (cellflag)
    [p, m] = size (sys);
    l = length (s);
    H = mat2cell (H, p, m, ones (1, l))(:);
  endif

endfunction


function H = __fracval__ (num, den, s)

  # if s has less frequencies than order of num or den, use standard way
  len_s = length (size (s));
  if (len_s == 3)
    n_freq = size (s,3);
  else
    n_freq = size (s,2);
  endif

  if (n_freq <= max (length (den), length (num)) - 1)

    H = polyval (num, s) ./ polyval (den, s);

  else

    H = __polyval__ (num, s) ./ __polyval__ (den, s);

  endif

endfunction


function p_val = __polyval__ (p, s)

  ## The code below of an alternative method for calculating
  ## a polynomial value does not seem to work in case of polynomial
  ## zeros in the origin. Therefore, separate the poylomials s^i
  ## from the original polynomial.
  if (p(end) == 0)
    p_red = flip (polyreduce (flip (p)));
    e = length (p) - length (p_red);
  else
    p_red = p;
    e = 0;
  endif

  ## Code suggested by 	dasergatskov in bug #63393, comment #17
  ## as numerically more robust variant for polyval, see
  ## https://savannah.gnu.org/bugs/index.php?63393
  p_val = polyval (p_red, s);
  [p_new, ~, mu] = polyfit (s, p_val, length (p_red)-1);
  p_val = polyval (p_new, s, [], mu) .* s.^e;

endfunction

%!test
%! num = [0, 0, 2.052877715426585e9, 4.416784109977014e13, 1.719831064785571e17];
%! den = [1, 8.380906856780281e4, 2.269817212624148e8, 4.344755797517798e12, 0];
%! H = freqresp (tf (num, den), [0, 100]);
%! assert (isinf (H(1, 1, 1)));
%! assert (imag (H(1, 1, 2)) < 0);
