## Copyright (C) 2009-2015   Lukas F. Reichlin
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
## Common code for adjusting TF model data.
## Used by tf and __set__.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2015
## Version: 0.1

function [num, den, tsam, tfvar] = __adjust_tf_data__ (num, den, tsam = -2)

  if (is_real_matrix (num) && isempty (den))    # static gain  tf (4),  tf (matrix)
    num = num2cell (num);
    den = num2cell (ones (size (num)));
    tsam = -2;
  endif

  if (! iscell (num))
    num = {num};
  endif
  
  if (! iscell (den))
    den = {den};
  endif

  num = cellfun (@tfpoly, num, "uniformoutput", false);
  den = cellfun (@tfpoly, den, "uniformoutput", false);
  
  if (tsam == 0)
    tfvar = "s";
  elseif (tsam == -2)
    tfvar = "x";
  else
    tfvar = "z";
  endif

endfunction
