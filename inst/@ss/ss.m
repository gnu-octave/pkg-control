## Copyright (C) 2009, 2010   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{sys} =} ss (@var{sys})
## @deftypefnx {Function File} {@var{sys} =} ss (@var{d})
## @deftypefnx {Function File} {@var{sys} =} ss (@var{a}, @var{b}, @var{c})
## @deftypefnx {Function File} {@var{sys} =} ss (@var{a}, @var{b}, @var{c}, @var{d}, @dots{})
## @deftypefnx {Function File} {@var{sys} =} ss (@var{a}, @var{b}, @var{c}, @var{d}, @var{tsam}, @dots{})
## Create or convert to state-space model.
##
## @strong{Inputs}
## @table @var
## @item a
## State transition matrix (n-by-n).
## @item b
## Input matrix (n-by-m).
## @item c
## Measurement matrix (p-by-n).
## @item d
## Feedthrough matrix (p-by-m).
## @item tsam
## Sampling time. If @var{tsam} is not specified, a continuous-time
## model is assumed.
## @item @dots{}
## Optional pairs of properties and values.
## Type @command{set (ss)} for more information.
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## State-space model.
## @end table
##
## @seealso{tf, dss}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.3

function sys = ss (a = [], b = [], c = [], d = [], varargin)

  ## model precedence: frd > ss > zpk > tf > double
  ## inferiorto ("frd");
  superiorto ("zpk", "tf", "double");

  argc = 0;                           # initialize argument count
  tsam = 0;                           # initialize sampling time

  switch (nargin)
    case {0, 3, 4}                    # ss (), ss (a, b, c), ss (a, b, c, d)
    ## nothing is done here
    ## case needed to prevent "otherwise"

    case 1
      if (isa (a, "ss"))              # already in ss form
        sys = a;
        return;
      elseif (isa (a, "lti"))         # another lti object
        [sys, alti] = __sys2ss__ (a);
        sys.lti = alti;               # preserve lti properties
        return;
      elseif (is_real_matrix (a))     # static gain  sys = ss (5)
        d = a;
        a = [];
      else
        print_usage ();
      endif

    case 2
      print_usage ();

    otherwise                         # default case  sys = ss (a, b, c, d, "prop1, "val1", ...)
      argc = numel (varargin);        # number of additional arguments after d
      if (issample (varargin{1}, 0))  # sys = ss (a, b, c, d, tsam, "prop1, "val1", ...)
        tsam = varargin{1};
        argc--;
        if (argc > 0)
          varargin = varargin(2:end);
        endif
      endif
  endswitch

  [a, b, c, d, tsam] = __adjust_ss_data__ (a, b, c, d, tsam);
  [p, m, n] = __ss_dim__ (a, b, c, d);

  stname = repmat ({""}, n, 1);

  ssdata = struct ("a", a, "b", b,
                   "c", c, "d", d,
                   "e", [],
                   "stname", {stname});

  ltisys = lti (p, m, tsam);

  sys = class (ssdata, "ss", ltisys);

  if (argc > 0)
    sys = set (sys, varargin{:});
  endif

endfunction
