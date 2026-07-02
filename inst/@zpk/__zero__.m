## Copyright (C) 2026  Mitchell Thompkins
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
## Return zeros of ZPK model directly from stored data (no polynomial round-trip).

## Author: Mitchell Thompkins <mitchell.thompkins@pm.me>
## Created: June 2026
## Version: 0.1

function [zer, gain, info] = __zero__ (sys, ~)

  if (issiso (sys))
    zer = sys.z{1};
    gain = sys.k;
    info = [];
  else
    warning ("Control:convert-to-state-space",
             "zpk: zero: converting to minimal state-space for zero computation of mimo zpk\n");
    [zer, gain, info] = zero (ss (sys));
  endif

endfunction
