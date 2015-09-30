## Copyright (C) 2009-2015   Lukas F. Reichlin
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
## @deftypefn{Function File} {@var{C} =} pidstd (@var{Kp})
## @deftypefnx{Function File} {@var{C} =} pidstd (@var{Kp}, @var{Ti})
## @deftypefnx{Function File} {@var{C} =} pidstd (@var{Kp}, @var{Ti}, @var{Td})
## @deftypefnx{Function File} {@var{C} =} pidstd (@var{Kp}, @var{Ti}, @var{Td}, @var{N})
## Return the transfer function @var{C} of the @acronym{PID} controller
## in standard form with first-order roll-off.
##
## @example
## @group
##                 1        Td s
## C(s) = Kp (1 + ---- + ----------)
##                Ti s   Td/N s + 1
## @end group
## @end example
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2015
## Version: 0.1

function C = pidstd (Kp = 1, Ti = inf, Td = 0, N = inf)

  if (! is_real_scalar (Kp, Ti, Td, N) || nargin > 4)
    print_usage ();
  endif

  if (N == inf)
    C = tf ([Kp*Td*Ti, Kp*Ti, Kp], [Ti, 0]);
  else
    C = tf ([Kp*N*Td*Ti+Kp*Td*Ti, Kp*N*Ti+Kp*Td, Kp*N], [Td*Ti, N*Ti, 0]);
  endif

endfunction
