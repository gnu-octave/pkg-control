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
## @deftypefn {Function File} {@var{retsys} =} xperm (@var{sys}, @var{idx})
## Reorder states in state-space models.
##
## @strong{Inputs}
## @table @var
## @item sys
## State-space model.
## @item idx
## Vector containing the state indices in the desired order.
## Alternatively, a cell vector containing the state names
## is possible as well.  See @code{sys.stname}.  State names
## only work if they were assigned explicitly before, i.e.
## @code{sys.stname} contains no empty strings.
## Note that if certain state indices of @var{sys} are
## missing or appear multiple times in @var{idx}, these
## states will be pruned or duplicated accordingly in the
## resulting state-space model @var{retsys}.
## @end table
##
## @strong{Outputs}
## @table @var
## @item retsys
## Resulting state-space model with states reordered according to @var{idx}.
## @end table
##
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: November 2009
## Version: 0.2

function sys = xperm (sys, st_idx)

  if (nargin != 2)
    print_usage ();
  endif

  if (! isa (sys, "ss"))
    warning ("xperm: system not in state-space form");
    sys = ss (sys);
  endif

  sys = __sys_prune__ (sys, ":", ":", st_idx);

endfunction
