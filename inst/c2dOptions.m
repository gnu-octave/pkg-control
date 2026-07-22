## Copyright (C) 2026        Prateek Ganguli
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
## @deftypefn {Function File} {@var{opt} =} c2dOptions (@var{'key1'}, @var{value1}, @var{'key2'}, @var{value2}, @dots{})
## Create c2d options struct @var{opt} from a number of key and value pairs.
## The struct contains configuration for continuous-to-discrete-time conversion.
##
## @strong{Properties}
## @table @var
## @item Method
## Discretization method (default: "zoh").
## @item PrewarpFrequency
## Prewarping frequency (default: 0).
## @item DelayModeling
## How to handle delays: "delay" or "state" (default: "delay").
## @end table
##
## @strong{Outputs}
## @table @var
## @item opt
## Struct with fields Method, PrewarpFrequency, DelayModeling.
## @end table
## @end deftypefn

## Author: Prateek Ganguli <prateek.ganguli@gmail.com>
## Created: July 2026
## Version: 0.1

function opt = c2dOptions (varargin)

  opt = struct ("Method", "zoh", "PrewarpFrequency", 0, "DelayModeling", "delay");

  if (rem (numel (varargin), 2))
    error ("c2dOptions: properties and values must come in pairs");
  endif

  for k = 1 : 2 : numel (varargin)
    key = varargin{k};
    val = varargin{k+1};

    switch (key)
      case "Method"
        if (! ischar (val))
          error ("c2dOptions: Method must be a string");
        endif
        opt.Method = val;

      case "PrewarpFrequency"
        if (! issample (val, 0))
          error ("c2dOptions: PrewarpFrequency must be a valid pre-warping frequency");
        endif
        opt.PrewarpFrequency = val;

      case "DelayModeling"
        val = lower (val);
        if (! any (strcmp (val, {"delay", "state"})))
          error ("c2dOptions: DelayModeling must be 'delay' or 'state'");
        endif
        opt.DelayModeling = val;

      otherwise
        error ("c2dOptions: unknown property '%s'", key);
    endswitch
  endfor

endfunction


%!test  # defaults
%! opt = c2dOptions ();
%! assert (opt.Method, "zoh");
%! assert (opt.PrewarpFrequency, 0);
%! assert (opt.DelayModeling, "delay");

%!test  # override DelayModeling only
%! opt = c2dOptions ("DelayModeling", "state");
%! assert (opt.DelayModeling, "state");
%! assert (opt.Method, "zoh");
%! assert (opt.PrewarpFrequency, 0);

%!test  # case-insensitive DelayModeling
%! opt = c2dOptions ("DelayModeling", "STATE");
%! assert (opt.DelayModeling, "state");

%!test  # override Method and PrewarpFrequency
%! opt = c2dOptions ("Method", "tustin", "PrewarpFrequency", 10);
%! assert (opt.Method, "tustin");
%! assert (opt.PrewarpFrequency, 10);

%!error <DelayModeling> c2dOptions ("DelayModeling", "bogus")
%!error <unknown property> c2dOptions ("NotAKey", 1)
%!error c2dOptions ("Method")
