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
## @deftypefn {Function File} {@var{db} =} mag2db (@var{mag})
## Convert Magnitude to Decibels (dB).
##
## @strong{Inputs}
## @table @var
## @item mag
## Magnitude value(s).  Both real-valued scalars and matrices are accepted.
## @end table
##
## @strong{Outputs}
## @table @var
## @item db
## Decibel (dB) value(s).
## @end table
##
## @seealso{db2mag}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: November 2012
## Version: 0.1

function db = mag2db (mag)

  if (nargin != 1 || ! is_real_matrix (mag))
    print_usage ();
  endif

  db = 20 .* log10 (mag);
  db(mag < 0) = NaN;

endfunction


%!assert (mag2db (100), 40);
%!assert (mag2db (0.1), -20);
