## Copyright (C) 2010   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{sys} =} frd (@var{sys})
## @deftypefnx {Function File} {@var{sys} =} frd (@var{sys}, @var{w})
## @deftypefnx {Function File} {@var{sys} =} frd (@var{H}, @var{w}, @dots{})
## @deftypefnx {Function File} {@var{sys} =} frd (@var{H}, @var{w}, @var{tsam}, @dots{})
## Create or convert to frequency response data.
##
## @strong{Inputs}
## @table @var
## @item sys
## LTI model.  If second argument @var{w} is omitted, the interesting
## frequency range is calculated by the zeros and poles of @var{sys}.
## @item H
## Frequency response array (p-by-m-by-lw).  In the SISO case,
## a vector (lw-by-1) or (1-by-lw) is accepted as well.
## @item w
## Frequency vector (lw-by-1) in radian per second [rad/s].
## @item tsam
## Sampling time.  If @var{tsam} is not specified, a continuous-time
## model is assumed.
## @item @dots{}
## Optional pairs of properties and values.
## Type @command{set (frd)} for more information.
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## Frequency response data object.
## @end table
##
## @seealso{dss, ss, tf}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: February 2010
## Version: 0.1

function sys = frd (H = [], w = [], varargin)

  ## NOTE: * There's no such thing as a static gain
  ##         because FRD objects are measurements,
  ##         not models.
  ##       * If something like  sys1 = frd (5)  existed,
  ##         it would cause troubles in cases like
  ##         sys2 = ss (...), sys = sys1 * sys2
  ##         because sys2 needs to be converted to FRD,
  ##         but sys1 contains no valid frequencies.
  ##       * Because of the reasons given above,
  ##         tsam = -2 and w = -1 don't make sense

  ## model precedence: frd > ss > zpk > tf > double
  superiorto ("ss", "zpk", "tf", "double");
  
  argc = 0;                             # initialize argument count

  switch (nargin)
    case 0                              # empty object  sys = frd ()
      tsam = 0;                         # tsam = -2  is *not* possible

    case 1
      if (isa (H, "frd"))               # already in frd form  sys = frd (frdsys)
        sys = H;
        return;                         # nothing more to do here
      elseif (isa (H, "lti"))           # another lti object  sys = frd (sys)
        [sys, alti] = __sys2frd__ (H);
        sys.lti = alti;                 # preserve lti properties
        return;                         # nothing more to do here
      else                              # sys = frd (H)  *must* fail
        print_usage ();
      endif

    case 2
      if (isa (H, "lti"))               # another lti object  sys = frd (sys, w)
        [sys, alti] = __sys2frd__ (H, w);
        sys.lti = alti;                 # preserve lti properties
        return;                         # nothing more to do here
      else                              # sys = frd (H, w)
        tsam = 0;                       # continuous-time
      endif

    otherwise                           # default case
      argc = numel (varargin);          # number of additional arguments after H and w
      if (issample (varargin{1}, -10))  # sys = frd (H, w, tsam, "prop1", val1, ...)
        tsam = varargin{1};             # discrete-time
        argc--;                         # tsam is not a property-value pair
        if (argc > 0)                   # if there are any properties and values ...
          varargin = varargin(2:end);   # remove tsam from property-value list
        endif
      else                              # sys = frd (H, w, "prop1", val1, ...)
        tsam = 0;                       # continuous-time
      endif
  endswitch

  [H, w] = __adjust_frd_data__ (H, w);
  [p, m] = __frd_dim__ (H, w);

  frdata = struct ("H", H, "w", w);     # struct for frd-specific data
  ltisys = lti (p, m, tsam);            # parent class for general lti data

  sys = class (frdata, "frd", ltisys);  # create frd object

  if (argc > 0)                         # if there are any properties and values, ...
    sys = set (sys, varargin{:});       # use the general set function
  endif

endfunction
