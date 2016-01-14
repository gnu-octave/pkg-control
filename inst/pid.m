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
## @deftypefn{Function File} {@var{C} =} pid (@var{Kp})
## @deftypefnx{Function File} {@var{C} =} pid (@var{Kp}, @var{Ki})
## @deftypefnx{Function File} {@var{C} =} pid (@var{Kp}, @var{Ki}, @var{Kd})
## @deftypefnx{Function File} {@var{C} =} pid (@var{Kp}, @var{Ki}, @var{Kd}, @var{Tf})
## Return the transfer function @var{C} of the @acronym{PID} controller
## in parallel form with first-order roll-off.
##
## @example
## @group
##              Ki      Kd s
## C(s) = Kp + ---- + --------
##              s     Tf s + 1
## @end group
## @end example
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: June 2015
## Version: 0.1

## TODO: discrete-time case, dozens of options, ...
##       If you wish to kill time with this repetitive task,
##       I'm happy to add your work :-)

function C = pid (Kp = 1, Ki = 0, Kd = 0, Tf = 0)

  if (! is_real_scalar (Kp, Ki, Kd, Tf) || nargin > 4)
    print_usage ();
  endif

  if (Kd == 0)    # catch cases like  pid (2, 0, 0, 3)
    Tf = 0;
  endif

  if (Ki == 0)    # minimal realization if  num(3) == 0  and  den(3) == 0
    C = tf ([Kp*Tf+Kd, Kp], [Tf, 1]);
  else
    C = tf ([Kp*Tf+Kd, Kp+Ki*Tf, Ki], [Tf, 1, 0]);
  endif

endfunction
