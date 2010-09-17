## Copyright (C) 2009 - 2010   Lukas F. Reichlin
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
## Minimal realization of SS models. The physical meaning of states is lost.
## Uses SLICOT TB01PD by courtesy of NICONET e.V. <http://www.slicot.org>

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.2

function retsys = __minreal__ (sys, tol)

  if (! isempty (sys.e))
    error ("ss: minreal: dss models not supported yet");
  endif

  if (tol == "def")
    tol = 0;
  elseif (tol > 1)
    error ("ss: minreal: require tol <= 1");
  endif

  [a, b, c, nr] = sltb01pd (sys.a, sys.b, sys.c, tol);

  retsys = ss (a, b, c, sys.d);
  retsys.lti = sys.lti;  # retain i/o names and tsam

endfunction
