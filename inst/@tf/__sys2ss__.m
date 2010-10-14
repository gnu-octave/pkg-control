## Copyright (C) 2009   Lukas F. Reichlin
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
## TF to SS conversion.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.1

function [retsys, retlti] = __sys2ss__ (sys)

  if (! issiso (sys))
    error ("tf: tf2ss: MIMO case not implemented yet");
  endif

  [num, den] = tfdata (sys);

  num = num{1, 1};
  den = den{1, 1};

  ## tfpoly ensures that there are no leading zeros
  if (length (num) > length (den))
    error ("tf: tf2ss: system must be proper");
  endif

  [a, b, c, d] = __tf2ss__ (num, den);

  retsys = ss (a, b, c, d);
  retlti = sys.lti;   # preserve lti properties

endfunction