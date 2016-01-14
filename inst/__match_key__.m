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
## Check whether string @var{str} matches one of the candidates
## in cell @var{props} exactly and return the matching property/key
## name.  If no exact match is found, check whether there is a
## candidate name in @var{props} which starts with @var{str} and
## return the partial match.  In case of ambiguity or a mismatch,
## raise an error.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: June 2015
## Version: 0.1

function key = __match_key__ (str, props, caller = "match_key")

  if (! ischar (str))
    error ("%s: key name must be a string", caller);
  endif

  ## exact matching - needed for e.g. iddata properties "u" and "userdata"
  idx = strcmpi (str, props);
  n = sum (idx);
  
  if (n == 1)           # 1 exact match
    key = lower (str);
    return;
  elseif (n > 1)        # props are not unique, this would be a bug in the control package
    error ("%s: key name '%s' is ambiguous", caller, str);
  endif
  
  ## partial matching - n was zero
  idx = strncmpi (str, props, length (str));
  n = sum (idx);
  
  if (n == 1)
    key = lower (props{idx});
    return;
  elseif (n > 1)
    error ("%s: key name '%s' is ambiguous", caller, str);
  endif
  
  error ("%s: key name '%s' is unknown", caller, str);
  
endfunction
