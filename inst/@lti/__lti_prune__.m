## Copyright (C) 2009   Lukas F. Reichlin
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
## Version: 0.1

function lti = __lti_prune__ (lti, out_idx, in_idx)

  m = numel (lti.inname);       # size before pruning!
  p = numel (lti.outname);

  lti.inname = lti.inname(in_idx);
  lti.outname = lti.outname(out_idx);
  
  lti.ingroup = structfun (@(x) sparse (x, 1:length(x), 1, m, length(x)), lti.ingroup, "uniformoutput", false);
  lti.outgroup = structfun (@(x) sparse (x, 1:length(x), 1, p, length(x)), lti.outgroup, "uniformoutput", false);

  lti.ingroup = structfun (@(x) x(in_idx, :), lti.ingroup, "uniformoutput", false);
  lti.outgroup = structfun (@(x) x(out_idx, :), lti.outgroup, "uniformoutput", false);
  
  [lti.ingroup, ~] = structfun (@find, lti.ingroup, "uniformoutput", false);
  [lti.outgroup, ~] = structfun (@find, lti.outgroup, "uniformoutput", false);

endfunction
