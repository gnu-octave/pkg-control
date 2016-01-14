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

## TODO: discrete-time case, dozens of options, ...
## NOTE: I don't see any need to implement 'pid' and 'pidstd'
##       as LTI classes like the 'dark side' does.
##       Returning a transfer function seems to be sufficient.
##       These functions' sole purpose is to help novice users
##       running their scripts with as little changes as possible.

function C = pidstd (Kp = 1, Ti = inf, Td = 0, N = inf)

  if (! is_real_scalar (Kp, Ti, Td, N) || nargin > 4)
    print_usage ();
  endif

  if (Kp == 0)                              # pidstd (0)
    C = tf (0);
  elseif (Ti == inf && N == inf)            # pidstd (Kp, inf, Td)
    C = tf (Kp*[Td, 1], [1]);
  elseif (Ti == inf)                        # pidstd (Kp, inf, Td, N)
    C = tf (Kp*[N*Td+Td, N], [Td, N]);
  elseif (N == inf)                         # pidstd (Kp, Ti, Td),  pidstd (Kp, Ti)
    C = tf (Kp*[Td*Ti, Ti, 1], [Ti, 0]);
  else                                      # pidstd (Kp, Ti, Td, N)
    C = tf (Kp*[N*Td*Ti+Td*Ti, N*Ti+Td, N], [Td*Ti, N*Ti, 0]);
  endif

endfunction
