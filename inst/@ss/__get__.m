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
## Access key values of SS objects.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.3

function val = __get__ (sys, key)

  switch (key)   # {<internal name>, <user name>}
    case "a"
      val = sys.a;
    case "b"
      val = sys.b;
    case "c"
      val = sys.c;
    case "d"
      val = sys.d;
    case "e"
      val = sys.e;
    case {"stname", "statename"}
      val = sys.stname;
    case "scaled"
      val = sys.scaled;
    otherwise
      error ("ss: get: invalid key name '%s'", key);
  endswitch

endfunction
