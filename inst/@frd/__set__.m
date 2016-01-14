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
## Set or modify keys of FRD objects.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2010
## Version: 0.2

function sys = __set__ (sys, key, val)

  switch (key)   # {<internal name>, <user name>}
    case {"h", "response"}
      val = __adjust_frd_data__ (val, sys.w);
      __frd_dim__ (val, sys.w);
      sys.H = val;
    case {"w", "frequency"}
      [~, val] = __adjust_frd_data__ (sys.H, val);
      __frd_dim__ (sys.H, val);
      sys.w = val;
    otherwise
      error ("frd: set: invalid key name '%s'", key);
  endswitch

endfunction
