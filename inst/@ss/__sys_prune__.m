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
## Submodel extraction and reordering for SS objects.
## This file is part of the Model Abstraction Layer.
## For internal use only.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.3

function sys = __sys_prune__ (sys, out_idx, in_idx, st_idx = ":")

  [sys.lti, out_idx, in_idx] = __lti_prune__ (sys.lti, out_idx, in_idx);
  
  if (ischar (st_idx) && ! strcmp (st_idx, ":"))
    st_idx = {st_idx};
  endif
  if (iscell (st_idx))
    st_idx = cellfun (@(x) __str2idx__ (sys.stname, x), st_idx);
  endif

  sys.a = sys.a(st_idx, st_idx);
  sys.b = sys.b(st_idx, in_idx);
  sys.c = sys.c(out_idx, st_idx);
  sys.d = sys.d(out_idx, in_idx);

  if (! isempty (sys.e))
    sys.e = sys.e(st_idx, st_idx);
  endif

  sys.stname = sys.stname(st_idx);

endfunction


## NOTE: This local '__str2idx__' function is different from
##       the one defined in the file '__str2idx__.m'.  Why?
##       There are input and output groups, but no state groups.
##       At least the 'dark side' does not have state groups.
##       However, I do contemplate 'stgroup'.
function idx = __str2idx__ (name, str)

  tmp = strcmp (name, str)(:);

  switch (nnz (tmp))
    case 1
      idx = find (tmp);
    case 0
      error ("ss: state name '%s' not found", str);
    otherwise
      error ("ss: state name '%s' is ambiguous", str);
  endswitch

endfunction
