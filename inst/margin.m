## Copyright (C) 2009 Doug Stewart and Lukas Reichlin
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

## -*- texinfo -*-
## @deftypefn{Function File} {[@var{gamma}, @var{phi}, @var{w_gamma}, @var{w_phi}] =} margin (@var{sys})
## Gain and phase margins.
##
## @strong{Inputs}
## @table @var
## @item sys
## Continuous time system.
## @end table
##
## @strong{Outputs}
## @table @var
## @item gamma
## Gain margin (as gain, not dbs).
## @item phi
## Phase margin (in degrees).
## @item w_gamma
## Radian frequency for the gain margin.
## @item w_phi
## Radian frequency for the phase margin.
## @end table
##
## @seealso{bode}
## @end deftypefn

## Version: 0.1alpha


function [gamma, phi, w_gamma, w_phi] = margin(sys)

  if (nargin != 1 || (! isstruct (sys)))
    print_usage ();
  endif

  ## Get Frequency Response
  [mag, pha, w] = bode (sys);

  ## fix the phase wrap around
  phw = unwrap (pha * pi/180);
  pha = phw * 180/pi;

  ## find the Gain Margin and its location
  [x, ix] = min (abs (pha + 180));

  if (x > 1)
    gamma = "INF";
    ix = length (pha);
  else
    gamma = 1 / mag(ix);
  endif

  w_gamma = w(ix);

  ## find the Phase Margin and its location
  [pmx, ipmx] = min (abs (mag - 1));

  phi = 180 + pha(ipmx);
  w_phi = w(ipmx);

endfunction 
