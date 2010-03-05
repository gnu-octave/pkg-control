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
## SS to TF conversion.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.1

function [retsys, retlti] = __sys2tf__ (sys)

  if (! issiso (sys))
    error ("ss: ss2tf: MIMO case not implemented yet");
  endif

  if (isempty (sys.a))  # static gain
    num = sys.d;
    den = 1;
  else  # default case
    [zer, gain] = zero (sys);
    
    num = gain * real (poly (zer));
    den = real (poly (sys.a));
  endif

  retsys = tf (num, den, get (sys, "tsam"));  # tsam needed to set appropriate tfvar
  retlti = sys.lti;   # preserve lti properties

endfunction