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
## @deftypefn {Function File} {} set (@var{sys})
## @deftypefnx {Function File} {} set (@var{sys}, @var{"key"}, @var{value}, @dots{})
## @deftypefnx {Function File} {@var{retsys} =} set (@var{sys}, @var{"key"}, @var{value}, @dots{})
## Set or modify properties of @acronym{LTI} objects.
## If no return argument @var{retsys} is specified, the modified @acronym{LTI} object is stored
## in input argument @var{sys}.  @command{set} can handle multiple properties in one call:
## @code{set (sys, 'key1', val1, 'key2', val2, 'key3', val3)}.
## @code{set (sys)} prints a list of the object's key names.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.7

function retsys = set (sys, varargin)

  if (nargin == 1)       # set (sys), sys = set (sys)

    [keys, vals] = __lti_keys__ (sys);
    nrows = numel (keys);

    str = strjust (strvcat (keys), "right");
    str = horzcat (repmat ("   ", nrows, 1), str, repmat (":  ", nrows, 1), strvcat (vals));

    disp (str);

    if (nargout != 0)    # function sys = set (sys, varargin)
      retsys = sys;      # would lead to unwanted output when using
    endif                # set (sys)

  else                   # set (sys, "key1", val1, ...), sys = set (sys, "key1", val1, ...)

    if (! isa (sys, "lti"));
      print_usage ();
    endif

    if (rem (nargin-1, 2))
      error ("lti: set: keys and values must come in pairs");
    endif

    [p, m] = size (sys);
    keys = __lti_keys__ (sys, true);

    for k = 1 : 2 : (nargin-1)
      key = __match_key__ (varargin{k}, keys, [class(sys), ": set"]);
      val = varargin{k+1};

      switch (key)
        case {"inname", "inputname"}
          sys.inname = __adjust_labels__ (val, m);

        case {"outname", "outputname"}
          sys.outname = __adjust_labels__ (val, p);

        case "tsam"
          if (issample (val, -1) && isdt (sys))
            sys.tsam = val;
          elseif (is_real_scalar (val) && val == 0 && isct (sys))
            sys.tsam = 0;
          elseif (is_real_matrix (val) && isempty (val) && isdt (sys))
            sys.tsam = -1;
          else
            error ("lti: set: invalid sampling time");
          endif

        case {"ingroup", "inputgroup"}
          if (isstruct (val) && all (size (val) == 1) ...
              && all (structfun (@(x) is_group_idx (x, m), val)))
            empty = structfun (@isempty, val);
            fields = fieldnames (val);
            sys.ingroup = rmfield (val, fields(empty));
          else
            error ("lti: set: key 'ingroup' requires a scalar struct containing valid input indices in the range [1, %d]", m);
          endif

        case {"outgroup", "outputgroup"}
          if (isstruct (val) && all (size (val) == 1) ...
              && all (structfun (@(x) is_group_idx (x, p), val)))
            empty = structfun (@isempty, val);
            fields = fieldnames (val);
            sys.outgroup = rmfield (val, fields(empty));
          else
            error ("lti: set: key 'outgroup' requires a scalar struct containing valid output indices in the range [1, %d]", p);
          endif

        case "name"
          if (ischar (val) && ndims (val) == 2 && (rows (val) == 1 || isempty (val)))
            sys.name = val;
          else
            error ("lti: set: key 'name' requires a string");
          endif

        case "notes"
          if (iscellstr (val))
            sys.notes = val;
          elseif (ischar (val))
            sys.notes = {val};
          else
            error ("lti: set: key 'notes' requires string or cell of strings");
          endif

        case "userdata"
          sys.userdata = val;

        case "lti"
          if (isa (val, "lti"))
            lti_keys = __lti_keys__ (val, false, false);
            n = numel (lti_keys);
            lti_vals = cell (n, 1);
            [lti_vals{1:n}] = get (val, lti_keys{:});
            for k = 1:n
              try
                sys = set (sys, lti_keys{k}, lti_vals{k});
              end_try_catch
            endfor
          else
            error ("lti: set: key 'lti' requires an LTI model");
          endif

        otherwise
          sys = __set__ (sys, key, val);
      endswitch
    endfor

    if (nargout == 0)    # set (sys, "key1", val1, ...)
      assignin ("caller", inputname (1), sys);
    else                 # sys = set (sys, "key1", val1, ...)
      retsys = sys;
    endif

  endif

endfunction


function bool = is_group_idx (idx, n)

  bool = (isempty (idx) || (is_real_vector (idx) && all (idx > 0) && all (idx <= n) && all (abs (fix (idx)) == idx)));

endfunction
