## Copyright (C) 2024   Torsten Lilge
##
## This file is part of GNU Octave's Control Package.
##
## This file free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This file is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with LTI Syncope.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## Remove leading zeros from a polynomial or a cell array of polynomials.
## polynomials of length one are not changed. For internal use only.

## Author: Torsten Lilge <ttl-octave@mailbox.org>

function p = __remove_leading_zeros__ (p)

  if (isa (p, "cell"))
    p = cellfun (@__remove_leading_zeros__, p, "uniformoutput", false);
    return;
  endif

  idx = find (p != 0);
  
  if (isempty (idx))
    return;
  else
    p = p(idx(1):end);
  endif

endfunction
