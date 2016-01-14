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
## @deftypefn {Function File} {} get (@var{sys})
## @deftypefnx {Function File} {@var{value} =} get (@var{sys}, @var{"key"})
## @deftypefnx {Function File} {[@var{val1}, @var{val2}, @dots{}] =} get (@var{sys}, @var{"key1"}, @var{"key2"}, @dots{})
## Access key values of @acronym{LTI} objects.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.3

function varargout = get (sys, varargin)

  if (nargin == 1)
    [keys, vals] = __lti_keys__ (sys);
    nrows = numel (keys);
    str = strjust (strvcat (keys), "right");
    str = horzcat (repmat ("   ", nrows, 1), str, repmat (":  ", nrows, 1), strvcat (vals));
    disp (str);
  else
    if (! isa (sys, "lti"))
      print_usage ();
    endif
    
    keys = __lti_keys__ (sys, true);
    
    for k = 1 : (nargin-1)
      key = __match_key__ (varargin{k}, keys, [class(sys), ": get"]);

      switch (key)
        case {"inname", "inputname"}
          val = sys.inname;
        case {"outname", "outputname"}
          val = sys.outname;
        case {"ingroup", "inputgroup"}
          val = sys.ingroup;
        case {"outgroup", "outputgroup"}
          val = sys.outgroup;
        case "tsam"
          val = sys.tsam;
        case "name"
          val = sys.name;
        case "notes"
          val = sys.notes;
        case "userdata"
          val = sys.userdata;
        otherwise
          val = __get__ (sys, key);
      endswitch

      varargout{k} = val;
    endfor
  endif

endfunction
