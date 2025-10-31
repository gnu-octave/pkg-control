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
    e = 1;
  endif

  ## Code suggested by 	dasergatskov in bug #63393, comment #17
  ## as numerically more robust variant for polyval, see
  ## https://savannah.gnu.org/bugs/index.php?63393
  p_val = polyval (p_red, s);
  [p_new, ~, mu] = polyfit (s, p_val, length (p_red)-1);
  p_val = polyval (p_new, s, [], mu) .* s.^e;

endfunction
