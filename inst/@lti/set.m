## Copyright (C) 2009 - 2010   Lukas F. Reichlin
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
## @deftypefn {Function File} set (@var{sys})
## @deftypefnx {Function File} {@var{sys} =} set (@var{sys}, @var{"property"}, @var{value})
## Set or modify properties of LTI objects.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.2

function retsys = set (sys, varargin)

  if (nargin == 1)  # set (sys), sys = set (sys)

    [propv, valv] = __propnames__ (sys);
    nrows = numel (propv);

    str = strjust (strvcat (propv), "right");
    str = horzcat (repmat ("   ", nrows, 1), str, repmat (":  ", nrows, 1), strvcat (valv));

    disp (str);

    if (nargout != 0)
      retsys = sys;
    endif

  elseif (nargout == 0)  # set (sys, "prop1", val1, ...)

    warning ("lti: get: use sys = get (sys, ""property1"", ...) to save changes");
    warning ("          octave does not support pass by reference");

  else  # sys = set (sys, "prop1", val1, ...)

    if (rem (nargin-1, 2))
      error ("lti: set: properties and values must come in pairs");
    endif

    [p, m] = size (sys);

    for k = 1 : 2 : (nargin-1)
      prop = lower (varargin{k});
      val = varargin{k+1};

      switch (prop)
        case {"inname", "inputname"}
          sys.inname = __checkname__ (val, m);

        case {"outname", "outputname"}
          sys.outname = __checkname__ (val, p);

        case {"tsam", "ts"}
          if (issample (val))
            sys.tsam = val;
            warning ("lti: set: use the editing of property ""%s"" with caution", prop);
            warning ("          it may lead to corrupted models");
          else
            error ("lti: set: invalid sampling time");
          endif

          ## TODO: use of c2d, d2c and d2d if tsam changes

        otherwise
          sys = __set__ (sys, prop, val);
      endswitch
    endfor

    retsys = sys;

  endif

endfunction