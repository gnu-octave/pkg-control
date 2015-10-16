## Copyright (C) 2009-2015   Lukas F. Reichlin
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
## Set or modify properties of TF objects.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.4

function sys = __set__ (sys, prop, val)

  switch (prop)  # {<internal name>, <user name>}
    case "num"
      num = __adjust_tf_data__ (val, sys.den);
      __tf_dim__ (num, sys.den);
      sys.num = num;

    case "den"
      [~, den] = __adjust_tf_data__ (sys.num, val);
      __tf_dim__ (sys.num, den);
      sys.den = den;

    case {"tfvar", "variable"}
      if (ischar (val))
        candidates = {"s", "p", "z", "q", "z^-1", "q^-1"};
        idx = strcmpi (val, candidates);
        if (any (idx))
          val = candidates{idx};
          if (isscalar (val))
            sys.tfvar = val;
            sys.inv = false;
          else
            sys.tfvar = val(1);
            sys.inv = true;
            ## FIXME: raise error if sys.inv and continuous-time
            ## FIXME: no switching between {"s", "p"} and {"z", "q", "z^-1", "q^-1"}
          endif
        else
          error ("tf: set: the string '%s' is not a valid transfer function variable", val);
        endif
      else
        error ("tf: set: property '%s' requires a string", prop);
      endif

    case "inv"
      if (! isdt (sys))
        error ("tf: set: property 'inv' requires discrete-time system");
      elseif (! isscalar (val))
        error ("tf: set: property 'inv' must be a scalar logical");
      else
        sys.inv = logical (val);
      endif

    otherwise
      error ("tf: set: invalid property name '%s'", prop);

  endswitch

endfunction
