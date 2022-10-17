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
## @deftypefn {Function File} {@var{sys} =} frd (@var{sys})
## @deftypefnx {Function File} {@var{sys} =} frd (@var{sys}, @var{w})
## @deftypefnx {Function File} {@var{sys} =} frd (@var{H}, @var{w}, @dots{})
## @deftypefnx {Function File} {@var{sys} =} frd (@var{H}, @var{w}, @var{tsam}, @dots{})
## Create or convert to frequency response data.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} model to be converted to frequency response data.
## If second argument @var{w} is omitted, the interesting
## frequency range is calculated by the zeros and poles of @var{sys}.
## @item H
## Frequency response array (p-by-m-by-lw).  H(i,j,k) contains the
## response from input j to output i at frequency k.  In the SISO case,
## a vector (lw-by-1) or (1-by-lw) is accepted as well.
## @item w
## Frequency vector (lw-by-1) in radian per second [rad/s].
## Frequencies must be in ascending order.
## @item tsam
## Sampling time in seconds.  If @var{tsam} is not specified,
## a continuous-time model is assumed.
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
## @strong{Option Keys and Values}
## @table @var
## @item 'H'
## Frequency response array.  See 'Inputs' for details.
##
## @item 'w'
## Frequency vector.  See 'Inputs' for details.
##
## @item 'tsam'
## Sampling time.  See 'Inputs' for details.
##
## @item 'inname'
## The name of the input channels in @var{sys}.
## Cell vector of length m containing strings.
## Default names are @code{@{'u1', 'u2', ...@}}
##
## @item 'outname'
## The name of the output channels in @var{sys}.
## Cell vector of length p containing strings.
## Default names are @code{@{'y1', 'y2', ...@}}
##
## @item 'ingroup'
## Struct with input group names as field names and
## vectors of input indices as field values.
## Default is an empty struct.
##
## @item 'outgroup'
## Struct with output group names as field names and
## vectors of output indices as field values.
## Default is an empty struct.
##
## @item 'name'
## String containing the name of the model.
##
## @item 'notes'
## String or cell of string containing comments.
##
## @item 'userdata'
## Any data type.
## @end table
##
## @seealso{dss, ss, tf}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: February 2010
## Version: 0.2

function sys = frd (varargin)

  ## NOTE: * There's no such thing as a static gain
  ##         because FRD objects are measurements,
  ##         not models.
  ##       * If something like  sys1 = frd (5)  existed,
  ##         it would cause troubles in cases like
  ##         sys2 = ss (...), sys = sys1 * sys2
  ##         because sys2 needs to be converted to FRD,
  ##         but sys1 contains no valid frequencies.
  ##       * However, things like  frd (ss (5))  should
  ##         be possible.

  ## model precedence: frd > ss > zpk > tf > double
  superiorto ("ss", "zpk", "tf", "double");

  if (nargin == 1 && isa (varargin{1}, "frd"))
    sys = varargin{1};
    return;
  elseif (nargin != 0 && nargin <= 2 ...
          && isa (varargin{1}, "lti"))
    [sys, lti] = __sys2frd__ (varargin{:});
    sys.lti = lti;                      # preserve lti properties
    return;
  endif

  H = []; w = [];                       # default frequency response data
  tsam = -2;                            # default sampling time
  
  [mat_idx, opt_idx, obj_flg] = __lti_input_idx__ (varargin);

  switch (numel (mat_idx))
    case 2                              # frd (H, w),  frd (H, w, ...)
      [H, w] = varargin{mat_idx};
      tsam = 0;
    case 3                              # frd (H, w, tsam),  frd (H, w, tsam, ...)
      [H, w, tsam] = varargin{mat_idx};
      if (isempty (tsam) && is_real_matrix (tsam))
        tsam = -1;
      elseif (! issample (tsam, -10))
        error ("frd: invalid sampling time");
      endif
    case 0                              # frd ()
      ## nothing to do here, just prevent case 'otherwise'
    otherwise                           # sys = frd (H)  *must* fail
      print_usage ();
  endswitch

  varargin = varargin(opt_idx);
  if (obj_flg)
    varargin = horzcat ({"lti"}, varargin);
  endif

  [H, w] = __adjust_frd_data__ (H, w);
  [p, m] = __frd_dim__ (H, w);          # determine number of outputs and inputs

  frdata = struct ("H", H, "w", w);     # struct for frd-specific data
  ltisys = lti (p, m, tsam);            # parent class for general lti data

  sys = class (frdata, "frd", ltisys);  # create frd object

  if (numel (varargin) > 0)             # if there are any properties and values, ...
    sys = set (sys, varargin{:});       # use the general set function
  endif

endfunction
