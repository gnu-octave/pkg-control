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
## Poles of SS object.

## Special thanks to Peter Benner for his advice.
## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.4

function pol = __pole__ (sys)

  if (isempty (sys.e))
    pol = eig (sys.a);
  else
    ## pol = eig (sys.a, sys.e);
    ## tol = norm ([sys.a, sys.e], 2);
    ## idx = find (abs (pol) < tol/eps);
    ## pol = pol(idx);
    ## do not scale, matrices B, C & D missing 
    pol = __sl_ag08bd__ (sys.a, sys.e, [], [], [], true);
  endif

endfunction
