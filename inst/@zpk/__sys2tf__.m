## Copyright (C) 2026  Mitchell Thompkins <mitchell.thompkins@pm.me>
##
## This file is part of the control package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not,
## see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## ZPK to TF conversion via poly().  Lossy for high-order systems; only used
## for explicit conversion (e.g. tf(zpk_sys)), not for c2d matched method.

## Author: Mitchell Thompkins <mitchell.thompkins@pm.me>
## Created: June 2026
## Version: 0.1

function [retsys, retlti] = __sys2tf__ (sys)

  [z, p, k] = __sys_data__ (sys);
  tsam = get (sys, "tsam");

  num = cellfun (@(zi, ki) real (ki * poly (zi)), z, num2cell (k), ...
                 "uniformoutput", false);
  den = cellfun (@(pi) real (poly (pi)), p, "uniformoutput", false);

  retsys = tf (num, den, tsam);
  retlti = sys.lti;

endfunction
