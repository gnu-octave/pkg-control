## Copyright (C) 2009, 2013   Lukas F. Reichlin
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
## Submodel extraction and reordering for LTI objects.
## This file is part of the Model Abstraction Layer.
## For internal use only.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.2

function lti = __lti_prune__ (lti, out_idx, in_idx)

  m = numel (lti.inname);       # size before pruning!
  p = numel (lti.outname);

  lti.inname = lti.inname(in_idx);
  lti.outname = lti.outname(out_idx);

  if (nfields (lti.ingroup) > 0)
    lti.ingroup = structfun (@(x) __group_prune__ (x, in_idx, m), lti.ingroup, "uniformoutput", false);
  endif
  if (nfields (lti.outgroup) > 0)
    lti.outgroup = structfun (@(x) __group_prune__ (x, out_idx, p), lti.outgroup, "uniformoutput", false);
  endif

endfunction


function group = __group_prune__ (group, idx, n)

  lg = length (group)
  group = sparse (group, 1:lg, 1, n, lg);
  group = group(idx, :);
  [group, ~] = find (group);

endfunction
