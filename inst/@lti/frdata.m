## Copyright (C) 2010   Lukas F. Reichlin
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
## @deftypefn {Function File} {[@var{H}, @var{w}, @var{tsam}] =} frdata (@var{sys})
## @deftypefnx {Function File} {[@var{H}, @var{w}, @var{tsam}] =} frdata (@var{sys}, @var{"vector"})
## Access frequency response data.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2010
## Version: 0.1

function [H, w, tsam] = frdata (sys, rtype = "vector")

  if (! isa (sys, "frd"))
    sys = frd (sys);
  endif

  [H, w] = __sys_data__ (sys); 

  tsam = sys.tsam;

  if (rtype(1) == "v")
    H = reshape (H, [], 1);
  endif

endfunction