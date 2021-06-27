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
## Access key values of TF objects.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.5

function val = __get__ (sys, key)

  switch (key)   # {<internal name>, <user name>}

    case {"num","den"}
      ## Give numerator and denominator of a tf component the same length
      [num, den] = __make_tf_polys_equally_long__ (sys);
      ## Get the coefficients of the polys with normalized length
      if key == "num"
        val = cellfun (@get, num, "uniformoutput", false);
      else
        val = cellfun (@get, den, "uniformoutput", false);
      endif

    case {"tfvar", "variable"}
      val = sys.tfvar;
      if (sys.inv && isdt (sys))
        val = [val, "^-1"];
      endif

    case "inv"
      val = sys.inv;

    otherwise
      error ("tf: get: invalid key name '%s'", key);

  endswitch

endfunction
