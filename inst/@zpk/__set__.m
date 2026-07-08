## Copyright (C) 2026  Mitchell Thompkins <mitchell.thompkins@pm.me>
##
## This file is part of the control package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not,
## see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## Set or modify keys of ZPK objects.

## Author: Mitchell Thompkins <mitchell.thompkins@pm.me>
## Created: June 2026
## Version: 0.1

function sys = __set__ (sys, key, val)

  switch (key)
    case {"z", "zeros"}
      if (! iscell (val))
        val = {val};
      endif
      if (! size_equal (val, sys.p, sys.k))
        error ("zpk: set: 'z' dimensions must match 'p' and 'k'");
      endif
      sys.z = cellfun (@(v) v(:), val, "uniformoutput", false);

    case {"p", "poles"}
      if (! iscell (val))
        val = {val};
      endif
      if (! size_equal (sys.z, val, sys.k))
        error ("zpk: set: 'p' dimensions must match 'z' and 'k'");
      endif
      sys.p = cellfun (@(v) v(:), val, "uniformoutput", false);

    case {"k", "gain"}
      if (! is_real_matrix (val))
        error ("zpk: set: 'k' must be a real-valued matrix");
      endif
      if (! size_equal (sys.z, sys.p, val))
        error ("zpk: set: 'k' dimensions must match 'z' and 'p'");
      endif
      sys.k = val;

    otherwise
      error ("zpk: set: invalid key name '%s'", key);
  endswitch

endfunction
