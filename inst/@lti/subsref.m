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
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## Subscripted reference for LTI objects.
## Used by Octave for "sys = sys(2:4, :)" or "val = sys.prop".

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.2

function a = subsref (a, s)

  if (isempty (s))
    error ("lti: subsref: missing index");
  endif

  for k = 1 : numel (s)
    if (isa (a, "lti"))
      switch (s(k).type)
        case "()"
          idx = s(k).subs;
          if (numel (idx) == 2)
            a = __sysprune__ (a, idx{1}, idx{2});
          elseif (numel (idx) == 1)
            a = __freqresp__ (a, idx{1});  
          else
            error ("lti: subsref: need one or two indices");
          endif
        case "."
          fld = s(k).subs;
          a = get (a, fld);
          ## warning ("lti: subsref: do not use subsref for development");
        otherwise
          error ("lti: subsref: invalid subscript type");
      endswitch
    else  # not an LTI model
      a = subsref (a, s(k:end));
      return;
    endif
  endfor

endfunction