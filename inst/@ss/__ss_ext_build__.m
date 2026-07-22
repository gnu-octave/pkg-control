## Copyright (C) 2026        Prateek Ganguli
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
## Build the extended ordinary system (A, [B1 B2], [C1;C2], [[D11 D12];[D21 D22]])
## used to discretize/continuize an InternalDelay-carrying ss model via the
## unchanged __c2d__/__d2c__ workers (which only look at matrix shapes).
##
## b2/c2/d12/d21/d22 are internal-only ss fields not exposed through
## get()/set(); this helper lives in @ss so it may access them directly
## (dot-access on an ss object from outside @ss is routed through the
## overloaded subsref/get, which only recognizes public keys).
##
## For internal use only, called from @lti/c2d.m and @lti/d2c.m.

## Author: Prateek Ganguli <prateek.ganguli@gmail.com>
## Created: July 2026
## Version: 0.1

function [ext_sys, nu, ny] = __ss_ext_build__ (sys)

  nu = columns (sys.b);
  ny = rows (sys.c);

  ext_sys = sys;
  ext_sys.b = [sys.b, sys.b2];
  ext_sys.c = [sys.c; sys.c2];
  ext_sys.d = [sys.d, sys.d12; sys.d21, sys.d22];

endfunction
