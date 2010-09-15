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
## Inversion of TF models.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.1

function sys = __sys_inverse__ (sys)

  [p, m] = size (sys);

  ## TODO: inversion of MIMO TF models

  if (p != 1 || m != 1)
    error ("tf: sys_inverse: MIMO systems not supported yet");
  endif

  den = sys.num;
  num = sys.den;

  if (den{1, 1} == 0)  # catch case den = 0
    num{1, 1} = tfpoly (0);
    den{1, 1} = tfpoly (1);
  endif

  sys.num = num;
  sys.den = den;

endfunction