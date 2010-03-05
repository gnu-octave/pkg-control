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
## @deftypefn {Function File} {@var{value} =} get (@var{sys}, @var{"property"})
## Access property values of LTI objects.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.1

function varargout = get (sys, varargin)

  if (nargin == 1)
    [propv, valv] = __propnames__ (sys);
    nrows = numel (propv);
    str = strjust (strvcat (propv), "right");
    str = horzcat (repmat ("   ", nrows, 1), str, repmat (":  ", nrows, 1), strvcat (valv));
    disp (str);
  else
    for k = 1 : (nargin-1)
      prop = lower (varargin{k});

      switch (prop)
        case {"inname", "inputname"}
          val = sys.inname;
        case {"outname", "outputname"}
          val = sys.outname;
        case {"tsam", "ts"}
          val = sys.tsam;
        otherwise
          val = __get__ (sys, prop);
      endswitch

      varargout{k} = val;
    endfor
  endif

endfunction