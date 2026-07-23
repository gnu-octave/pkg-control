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

## Author: Prateek Ganguli <prateek.ganguli@gmail.com>
## Created: July 2026
## Version: 0.1

## Delay-buffer lookup: w(port) = z_hist(k - tau(port), port), read as 0 when
## k - tau(port) < 1 (no past sample yet -- the zero-history boundary).  At
## k == tau this still reads 0 (index 0); the first real-history read is at
## k == tau + 1, returning z_hist(1) -- i.e. exactly a tau-sample delay,
## matching the integer-sample shift used by __apply_timeresp_delay__.
## Shared by __time_response__.m and lsim.m.
##
## Assumes tau(port) >= 1 for every port: a port whose delay rounds to 0
## samples (an internal delay smaller than half a sample time) would read
## w(k) = z_hist(k), which is always 0 at read time since z_hist(k) is only
## written later in the same step -- silently dropping that port's D12/D22
## feedthrough instead of solving the resulting algebraic w=z loop. Not
## reachable via a well-formed nonzero delay; c2d's rounding is expected to
## keep tau >= 1 for any InternalDelay this function is asked to simulate.
function w = __delay_lookup__ (z_hist, k, tau, nports)
  w = zeros (nports, 1);
  for pp = 1 : nports
    idx = k - tau(pp);
    if (idx >= 1)
      w(pp) = z_hist(idx, pp);
    endif
  endfor
endfunction
