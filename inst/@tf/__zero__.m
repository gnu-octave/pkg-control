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
## Transmission zeros of TF object.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.3

function [zer, gain, info] = __zero__ (sys)

  if (issiso (sys))
    num = get (sys.num{1});
    den = get (sys.den{1});
    zer = roots (num);
    gain = num(1) / den(1);
    info = [];
  else
    warning ("tf: zero: converting to minimal state-space for zero computation of mimo tf");
    [zer, gain, info] = zero (ss (sys));
  endif

endfunction
