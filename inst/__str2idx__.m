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
## Get input/output indices from in/outgroup and in/outname.
## For internal use only.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2013
## Version: 0.1

function idx = __str2idx__ (group, name, str, id)

  if (isfield (group, str))
    idx = group.(str)(:);
  else
    tmp = strcmp (name, str)(:);
    switch (nnz (tmp))
      case 1
        idx = find (tmp);
      case 0
        error ("lti: %sgroup or %sname '%s' not found", id, id, str);
      otherwise
        error ("lti: %sname '%s' is ambiguous", id, str);
        ## FIXME: error for structure arrays
    endswitch
  endif

endfunction
