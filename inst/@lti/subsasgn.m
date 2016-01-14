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
## Subscripted assignment for @acronym{LTI} objects.
## Used by Octave for "sys.property = value".

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.3

function sys = subsasgn (sys, idx, val)

  ## TODO: enable *more* stuff like sys.a(2, 1:3) = [4, 5, 6]
  ## warning ("lti: subsasgn: do not use subsasgn for development");

  switch (idx(1).type)
    case "."
      if (length (idx) == 1)
        sys = set (sys, idx.subs, val);
      else
        key = idx(1).subs;
        sys = set (sys, key, subsasgn (get (sys, key), idx(2:end), val));
      endif

    case "()"
      if (numel (idx(1).subs) > 2)      # sys(:, :, ...) = ltisys
        ## bug #45314  <https://savannah.gnu.org/bugs/?45314>
        error (["lti: subsasgn: arrays of linear models are not supported (yet).  ", ...
                "as a workaround, try to use cells of LTI models instead."]);
      else                              # sys(out_idx, in_idx).a = mat  
        error (["lti: subsasgn: invalid subscripted assignment type '()'.  ", ...
                "you must select an LTI key first, i.e.  sys.keyname(...) = ..."]);
      endif

    otherwise
      error ("lti: subsasgn: invalid subscripted assignment type '%s'", idx(1).type);

  endswitch

endfunction
