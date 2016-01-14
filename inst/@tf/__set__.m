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
## Set or modify keys of TF objects.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.4

function sys = __set__ (sys, key, val)

  switch (key)   # {<internal name>, <user name>}
    case "num"
      if (get (sys, "tsam") == -2)   # NOTE: sys.lti.tsam  would call 'get' via 'subsref'
        error (["tf: set: tinkering with numerators of static gains is disabled on purpose.  " ...
                "to avoid this error, set the sampling time of your LTI model first."]);
      elseif (sys.inv)               # FIXME: f = filt (1, [1 0 0 1]),  f.num = 1
        warning ("tf:num-inv", ...
                 "tf: although system is displayed in '%s^-1', numerator 'num' must be specified in '%s'", ...
                 sys.tfvar, sys.tfvar);
      endif
      num = __adjust_tf_data__ (val, sys.den);
      __tf_dim__ (num, sys.den);
      sys.num = num;

    case "den"
      if (get (sys, "tsam") == -2)
        error (["tf: set: tinkering with denominators of static gains is disabled on purpose.  " ...
                "to avoid this error, set the sampling time of your LTI model first."]);
      elseif (sys.inv)               # FIXME: f = filt ([1 0 0 1], 1),  f.den = 1
        warning ("tf:den-inv", ...
                 "tf: although system is displayed in '%s^-1', denominator 'den' must be specified in '%s'", ...
                 sys.tfvar, sys.tfvar);      
      endif
      [~, den] = __adjust_tf_data__ (sys.num, val);
      __tf_dim__ (sys.num, den);
      sys.den = den;

    case {"tfvar", "variable"}
      if (ischar (val))
        candidates = {"s", "p", "z", "q", "z^-1", "q^-1"};
        idx = strcmpi (val, candidates);
        if (any (idx))
          val = candidates{idx};
          n = find (idx);
          if (n > 2 && isct (sys))
            error ("tf: set: variable '%s' not allowed for static gains and continuous-time models", val);
          elseif (n < 3 && isdt (sys))
            error ("tf: set: variable '%s' not allowed for static gains and discrete-time models", val);
          endif
          if (isscalar (val))
            sys.tfvar = val;
            sys.inv = false;
          else
            sys.tfvar = val(1);
            sys.inv = true;
          endif
        else
          error ("tf: set: the string '%s' is not a valid transfer function variable", val);
        endif
      else
        error ("tf: set: key '%s' requires a string", prop);
      endif

    case "inv"
      if (islogical (val) && isscalar (val) || is_real_scalar (val))
        sys.inv = logical (val);
      else
        error ("tf: set: key 'inv' must be a scalar logical");
      endif

    otherwise
      error ("tf: set: invalid key name '%s'", key);
  endswitch

endfunction
