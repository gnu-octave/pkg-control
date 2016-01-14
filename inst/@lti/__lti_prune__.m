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
## Submodel extraction and reordering for LTI objects.
## This file is part of the Model Abstraction Layer.
## For internal use only.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.2

function [lti, out_idx, in_idx] = __lti_prune__ (lti, out_idx, in_idx)

  if (ischar (out_idx) && ! strcmp (out_idx, ":"))  # sys("grp", :)
    out_idx = {out_idx};
  endif
  if (ischar (in_idx) && ! strcmp (in_idx, ":"))    # sys(:, "grp")
    in_idx = {in_idx};
  endif

  if (iscell (out_idx))                             # sys({"grp1", "grp2"}, :)
    tmp = cellfun (@(x) __str2idx__ (lti.outgroup, lti.outname, x, "out"), out_idx, "uniformoutput", false);
    out_idx = vertcat (tmp{:});
  endif
  if (iscell (in_idx))                              # sys(:, {"grp1", "grp2"})
    tmp = cellfun (@(x) __str2idx__ (lti.ingroup, lti.inname, x, "in"), in_idx, "uniformoutput", false);
    in_idx = vertcat (tmp{:});
  endif

  if (nfields2 (lti.outgroup))
    p = numel (lti.outname);                        # get size before pruning outnames!
    [lti.outgroup, empty] = structfun (@(x) __group_prune__ (x, out_idx, p), lti.outgroup, "uniformoutput", false);
    empty = cell2mat (struct2cell (empty));
    fields = fieldnames (lti.outgroup);
    lti.outgroup = rmfield (lti.outgroup, fields(empty));
  endif
  if (nfields2 (lti.ingroup))
    m = numel (lti.inname); 
    [lti.ingroup, empty] = structfun (@(x) __group_prune__ (x, in_idx, m), lti.ingroup, "uniformoutput", false);
    empty = cell2mat (struct2cell (empty));
    fields = fieldnames (lti.ingroup);
    lti.ingroup = rmfield (lti.ingroup, fields(empty));
  endif
  
  lti.outname = lti.outname(out_idx);
  lti.inname = lti.inname(in_idx);

endfunction


function [group, empty] = __group_prune__ (group, idx, n)

  lg = length (group);
  group = sparse (group, 1:lg, 1, n, lg);
  group = group(idx, :);
  [group, ~] = find (group);
  empty = isempty (group);

endfunction
