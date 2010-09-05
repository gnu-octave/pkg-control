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
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## Number of inputs (m), states (n) and outputs (p) of state space matrices.
## For internal use only.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.2

function [m, n, p] = __ssmatdim__ (a, b, c, d)

  [arows, acols] = size (a);
  [brows, bcols] = size (b);
  [crows, ccols] = size (c);
  [drows, dcols] = size (d);

  m = bcols;  # = dcols
  n = arows;  # = acols
  p = crows;  # = drows

  if (! issquare (a) && ! isempty (a))
    error ("ss: system matrix a(%dx%d) is not square", arows, acols);
  endif

  if (brows != arows)
    error ("ss: system matrices a(%dx%d) and b(%dx%d) are incompatible",
            arows, acols, brows, bcols);
  endif

  if (ccols != acols)
    error ("ss: system matrices a(%dx%d) and c(%dx%d) are incompatible",
            arows, acols, crows, ccols);
  endif

  if (bcols != dcols)
    error ("ss: system matrices b(%dx%d) and d(%dx%d) are incompatible",
            brows, bcols, drows, dcols);
  endif

  if (crows != drows)
    error ("ss: system matrices c(%dx%d) and d(%dx%d) are incompatible",
            crows, ccols, drows, dcols);
  endif

  if (! isreal (a) || ! isreal (b) || ! isreal (c) || ! isreal (d))
    error ("ss: system matrices are not real");
  endif

endfunction
