## Copyright (C) 2009-2015   Lukas F. Reichlin
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
## @deftypefn {Function File} {} set (@var{dat})
## @deftypefnx {Function File} {} set (@var{dat}, @var{"property"}, @var{value}, @dots{})
## @deftypefnx {Function File} {@var{dat} =} set (@var{dat}, @var{"property"}, @var{value}, @dots{})
## Set or modify properties of iddata objects.
## If no return argument @var{dat} is specified, the modified @acronym{LTI} object is stored
## in input argument @var{dat}.  @command{set} can handle multiple properties in one call:
## @code{set (dat, 'prop1', val1, 'prop2', val2, 'prop3', val3)}.
## @code{set (dat)} prints a list of the object's property names.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: February 2012
## Version: 0.4

function retdat = set (dat, varargin)

  if (nargin == 1)       # set (dat), dat = set (dat)

    [props, vals] = __property_names__ (dat);
    nrows = numel (props);

    str = strjust (strvcat (props), "right");
    str = horzcat (repmat ("   ", nrows, 1), str, repmat (":  ", nrows, 1), strvcat (vals));

    disp (str);

    if (nargout != 0)       # function dat = set (dat, varargin)
      retdat = dat;         # would lead to unwanted output when using
    endif                   # set (dat)

  else                      # set (dat, "prop1", val1, ...), dat = set (dat, "prop1", val1, ...)

    if (! isa (dat, "iddata"))
      print_usage ();
    endif

    if (rem (nargin-1, 2))
      error ("iddata: set: properties and values must come in pairs");
    endif

    [n, p, m, e] = size (dat);
    keys = __property_names__ (dat, true);

    for k = 1 : 2 : (nargin-1)
      prop = __match_key__ (varargin{k}, keys, "iddata: set");
      val = varargin{k+1};

      switch (prop)
        case {"y", "outdata", "outputdata"}
          val = __adjust_iddata__ (val, dat.u);
          [pval, ~, eval] = __iddata_dim__ (val, dat.u);
          if (pval != p)
            error ("iddata: set: argument has %d instead of %d outputs", pval, p);
          endif
          if (eval != e)    # iddata_dim is not sufficient if dat.u = []
            error ("iddata: set: argument has %d instead of %d experiments", eval, e);
          endif
          if (dat.timedomain && ! is_real_matrix (val{:}))
            error ("iddata: set: require real-valued output signals for time domain datasets");
          endif
          dat.y = val;
        case {"u", "indata", "inputdata"}
          [~, val] = __adjust_iddata__ (dat.y, val);
          [~, mval] = __iddata_dim__ (dat.y, val);
          if (mval != m)
            error ("iddata: set: argument has %d instead of %d inputs", mval, m);
          endif
          if (dat.timedomain && ! is_real_matrix (val{:}))
            error ("iddata: set: require real-valued input signals for time domain datasets");
          endif
          dat.u = val;
        case {"outname", "outputname"}
          dat.outname = __adjust_labels__ (val, p);
        case {"inname", "inputname"}
          dat.inname = __adjust_labels__ (val, m);
        case {"outunit", "outputunit"}
          dat.outunit = __adjust_labels__ (val, p);
        case {"inunit", "inputunit"}
          dat.inunit = __adjust_labels__ (val, m);
        case {"timeunit"}
          if (ischar (val))
            dat.timeunit = val;
          else
            error ("iddata: set: property 'timeunit' requires a string");
          endif
        case {"expname", "experimentname"}
          dat.expname = __adjust_labels__ (val, e);
        case {"tsam"}
          dat.tsam = __adjust_iddata_tsam__ (val, e);
        case {"w", "frequency"}
          if (! iscell (val))
            val = {val};
          endif
          
          if (any (cellfun (@(w) ! isempty (w) && (! is_real_vector (w) || any (w < 0) ...
                                                   || ! issorted (w) || w(1) > w(end) ...
                                                   || length (unique (w)) != length (w)), val)))
            error ("iddata: set: w must be a vector of positive real values in ascending order");
          endif
          dat.w = val;
          dat.timedomain = false;
        case "name"
          if (ischar (val))
            dat.name = val;
          else
            error ("iddata: set: property 'name' requires a string");
          endif
        case "notes"
          if (iscellstr (val))
            dat.notes = val;
          elseif (ischar (val))
            dat.notes = {val};
          else
            error ("lti: set: property 'notes' requires string or cell of strings");
          endif
        case "userdata"
          dat.userdata = val;
        otherwise
          error ("iddata: set: invalid property name '%s'", varargin{k});
      endswitch
    endfor

    if (nargout == 0)    # set (dat, "prop1", val1, ...)
      assignin ("caller", inputname (1), dat);
    else                 # dat = set (dat, "prop1", val1, ...)
      retdat = dat;
    endif

  endif

endfunction