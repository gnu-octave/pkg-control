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
## Inversion of TF models.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.4

function sys = __sys_inverse__ (sys)

  if (issiso (sys))         # SISO
    num = sys.num;
    den = sys.den;
    if (num{1,1} == 0)      # catch case num = 0
      sys.num(1,1) = tfpoly (0);
      sys.den(1,1) = tfpoly (1);
    else
      sys.num = den;
      sys.den = num;
    endif
  else                      # MIMO
    ## I've calculated TF inversion of 2x2 and 3x3 systems with Sage CAS,
    ## but the formulae give systems with very high orders, therefore
    ## I always use the conversion to state-space and back.
    [num, den] = tfdata (inv (ss (sys)), "tfpoly");
    sys.num = num;
    sys.den = den;
  endif

endfunction
