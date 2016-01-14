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
## End indexing for @acronym{IDDATA} objects.
## Used by Octave for "dat(1:end)".

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: December 2013
## Version: 0.1

function ret = end (dat, k, n)

  if (n > 4)
    error ("iddata: end: require at most 4 indices in the expression");
  endif

  switch (k)
    case 1          # selecting samples
      ret = size (dat, 1);
      if (numel (ret) != 1 && ! isequal (num2cell (ret){:}))
        error ("iddata: end: for multi-experiment datasets, require equal number of samples when selecting samples with 'end'");
      endif
      ret = ret(1);
    case {2, 3, 4}  # selecting outputs, inputs or experiments
      ret = size (dat, k);
    otherwise
      error ("iddata: end: invalid expression index k = %d", k);
  endswitch

endfunction
