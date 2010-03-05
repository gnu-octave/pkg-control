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
## Version: 0.1

function ret = subsref (sys, s)

  if (isempty (s))
    error ("lti: subsref: missing index");
  endif

  switch (s(1).type)
    case "()"
      idx = s(1).subs;
      if (numel (idx) == 2)
        ret = __sysprune__ (sys, idx{1}, idx{2});
      elseif (numel (idx) == 1)
        ret = __freqresp__ (sys, idx{1});  
      else
        error ("lti: subsref: need one or two indices");
      endif

    case "."
      fld = s.subs;
      ret = get (sys, fld);
      ## warning ("lti: subsref: do not use subsref for development");

    otherwise
      error ("lti: subsref: invalid subscript type");

  endswitch

endfunction