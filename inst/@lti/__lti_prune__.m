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

function [lti, out_idx, in_idx] = __lti_prune__ (lti, out_idx, in_idx)

  if (ischar (out_idx) && ! strncmpi (out_idx, ":", 1))     # sys("grp", :)
    out_idx = lti.outgroup.(out_idx);
  elseif (iscell (out_idx))                                 # sys({"grp1", "grp2"}, :)
    tmp = cellfun (@(x) lti.outgroup.(x)(:), out_idx, "uniformoutput", false);
    out_idx = vertcat (tmp{:});
  endif

  if (ischar (in_idx) && ! strncmpi (in_idx, ":", 1))       # sys(:, "grp")
    in_idx = lti.ingroup.(in_idx);
  elseif (iscell (in_idx))                                  # sys(:, {"grp1", "grp2"})
    tmp = cellfun (@(x) lti.ingroup.(x)(:), in_idx, "uniformoutput", false);
    in_idx = vertcat (tmp{:});
  endif
  
  ## TODO: check group with isfield whether idx is field,
  ##       if not, also check inname/outname

  if (nfields (lti.outgroup))
    p = numel (lti.outname);                                # get size before pruning outnames!
    lti.outgroup = structfun (@(x) __group_prune__ (x, out_idx, p), lti.outgroup, "uniformoutput", false);
  endif
  if (nfields (lti.ingroup))
    m = numel (lti.inname); 
    lti.ingroup = structfun (@(x) __group_prune__ (x, in_idx, m), lti.ingroup, "uniformoutput", false);
  endif
  
  lti.outname = lti.outname(out_idx);
  lti.inname = lti.inname(in_idx);

endfunction


function group = __group_prune__ (group, idx, n)

  lg = length (group);
  group = sparse (group, 1:lg, 1, n, lg);
  group = group(idx, :);
  [group, ~] = find (group);

endfunction
