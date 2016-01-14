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
## Check whether a cell contains the required number of strings.
## Used by set and __set__.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.4

function name = __adjust_labels__ (name, req_len)

  if (iscell (name))
    name = reshape (name, [], 1);
  else                             # catch the siso case,
    name = {name};                 # e.g. sys = set (sys, "inname", "u_1")
  endif

  if (! iscellstr (name))
    error ("lti: set: require string or cell of strings");
  endif

  if (numel (name) != req_len)
    if (numel (name) == 1 && req_len > 1)
      if (isempty (name{1}))       # delete names quickly
        name = repmat ({""}, req_len, 1);
      else
        name = strseq (name{1}, 1:req_len);
      endif
    else
      error ("lti: set: cell must contain %d strings", req_len);
    endif
  endif

endfunction
