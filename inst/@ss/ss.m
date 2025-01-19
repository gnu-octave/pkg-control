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
## @deftypefn {Function File} {@var{sys} =} ss (@var{sys})
## @deftypefnx {Function File} {@var{sys} =} ss (@var{d}, @dots{})
## @deftypefnx {Function File} {@var{sys} =} ss (@var{a}, @var{b}, @dots{})
## @deftypefnx {Function File} {@var{sys} =} ss (@var{a}, @var{b}, @var{c}, @dots{})
## @deftypefnx {Function File} {@var{sys} =} ss (@var{a}, @var{b}, @var{c}, @var{d}, @dots{})
## @deftypefnx {Function File} {@var{sys} =} ss (@var{a}, @var{b}, @var{c}, @var{d}, @var{tsam}, @dots{})
## Create or convert to state-space model.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} model to be converted to state-space.
## @item a
## State matrix (n-by-n).
## @item b
## Input matrix (n-by-m).
## @item c
## Output matrix (p-by-n).
## If @var{c} is empty @code{[]} or not specified, an identity matrix is assumed.
## @item d
## Feedthrough matrix (p-by-m).
## If @var{d} is empty @code{[]} or not specified, a zero matrix is assumed.
## @item tsam
## Sampling time in seconds.  If @var{tsam} is not specified, a continuous-time model is assumed.
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
## @strong{Option Keys and Values}
## @table @var
## @item 'a', 'b', 'c', 'd', 'e'
## State-space matrices.  See 'Inputs' for details.
##
## @item 'stname'
## The name of the states in @var{sys}.
## Cell vector containing strings for each state.
## Default names are @code{@{'x1', 'x2', ...@}}
##
## @item 'scaled'
## Logical.  If set to true, no automatic scaling is used,
## e.g. for frequency response plots.
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
## @strong{Equations}
## @tex
## $$\dot{x} = A\,x + B\,u$$
## $$y = C\,x + D\,u$$
## @end tex
## @ifnottex
## @example
##   .
##   x = A x + B u
##   y = C x + D u
## @end example
## @end ifnottex
##
## @strong{Example}
## @example
## @group
## octave:1> a = [1 2 3; 4 5 6; 7 8 9];
## octave:2> b = [10; 11; 12];
## octave:3> stname = @{'V', 'A', 'kJ'@};
## octave:4> sys = ss (a, b, 'stname', stname)
##
## sys.a =
##         V   A  kJ
##    V    1   2   3
##    A    4   5   6
##    kJ   7   8   9
##
## sys.b =
##        u1
##    V   10
##    A   11
##    kJ  12
##
## sys.c =
##         V   A  kJ
##    y1   1   0   0
##    y2   0   1   0
##    y3   0   0   1
##
## sys.d =
##        u1
##    y1   0
##    y2   0
##    y3   0
##
## Continuous-time model.
## octave:5>
## @end group
## @end example
##
## @strong{Note on compatibility}
##
## If the state-space model @var{sys} is converted from a transfer
## function, the resulting state-space model can be transformed into
## the form computed by Matlab (a controllable canonical form with
## flipped state variables order) by using the following
## similarity transformation:
##
## @example
## @group
## n = size (sys.a, 1)
## QSi = inv (ctrb (sys))
## T(n,:) = QSi(n,:)
## for i=n-1:-1:1, T(i,:) = T(i+1,:)*sys.a, endfor
## sys_ml = ss2ss (sys, T)
## @end group
## @end example
##
## @seealso{tf, dss}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.4

function sys = ss (varargin)

  ## model precedence: frd > ss > zpk > tf > double
  ## inferiorto ("frd");
  superiorto ("zpk", "tf", "double");

  if (nargin == 1)                      # shortcut for lti objects
    if (isa (varargin{1}, "ss"))        # already in ss form  sys = ss (sssys)
      sys = varargin{1};
      return;
    elseif (isa (varargin{1}, "lti"))   # another lti object  sys = ss (sys)
      [sys, lti] = __sys2ss__ (varargin{1});
      sys.lti = lti;                    # preserve lti properties
      return;
    endif
  elseif (nargin == 2)                  # shortcut for lti objects plus
    if (ischar (varargin{2}) && isa (varargin{1}, "lti"))
      candidates = {"minimal", "explicit"};
      key = __match_key__ (varargin{2}, candidates, "ss");
      switch (key)
        case "minimal"                  # ss (sys, 'minimal')
          sys = minreal (ss (varargin{1}));
          return;
        case "explicit"                 # ss (sys, 'explicit')
          sys = dss2ss (ss (varargin{1}));
          return;
        otherwise                       # this would be a silly bug
      endswitch
    endif
  endif

  a = []; b = []; c = []; d = [];       # default state-space matrices
  tsam = 0;                             # default sampling time

  [mat_idx, opt_idx, obj_flg] = __lti_input_idx__ (varargin);

  switch (numel (mat_idx))
    case 1
      d = varargin{mat_idx};
    case 2
      [a, b] = varargin{mat_idx};
    case 3
      [a, b, c] = varargin{mat_idx};
    case 4
      [a, b, c, d] = varargin{mat_idx};
    case 5
      [a, b, c, d, tsam] = varargin{mat_idx};
      if (isempty (tsam) && is_real_matrix (tsam))
        tsam = -1;
      elseif (! issample (tsam, -10))
        error ("ss: invalid sampling time");
      endif
    case 0
      ## nothing to do here, just prevent case 'otherwise'
    otherwise
      print_usage ();
  endswitch

  varargin = varargin(opt_idx);
  if (obj_flg)
    varargin = horzcat ({"lti"}, varargin);
  endif

  [a, b, c, d, tsam] = __adjust_ss_data__ (a, b, c, d, tsam);
  [p, m, n] = __ss_dim__ (a, b, c, d);  # determine number of outputs, inputs and states

  stname = repmat ({""}, n, 1);         # cell with empty state names

  ssdata = struct ("a", a, "b", b,
                   "c", c, "d", d,
                   "e", [],
                   "stname", {stname},
                   "scaled", false);    # struct for ss-specific data

  ltisys = lti (p, m, tsam);            # parent class for general lti data

  sys = class (ssdata, "ss", ltisys);   # create ss object

  if (numel (varargin) > 0)             # if there are any properties and values, ...
    sys = set (sys, varargin{:});       # use the general set function
  endif

endfunction


## TODO: create a separate function @lti/dss2ss.m
function G = dss2ss (G)

  [G.a, G.b, G.c, G.d, G.e] = __dss2ss__ (G.a, G.b, G.c, G.d, G.e);

endfunction

%!test
%! # Test the example for the SLICOT function TD04AD
%! den1x = [1.0   6.0  11.0   6.0];
%! den2x = [1.0   6.0  11.0   6.0];
%! num11 = [1.0   6.0  12.0   7.0];
%! num12 = [0.0   1.0   4.0   3.0];
%! num21 = [0.0   0.0   1.0   1.0];
%! num22 = [1.0   8.0  20.0  15.0];
%! G = tf ( {num11, num12 ; num21, num22} , {den1x, den1x ; den2x, den2x} );
%! [A,B,C,D] = ssdata (ss (G));
%!
%! Ae = [ 0.499999999999998  -0.802831617715081   0.938714854383612 ;
%!        4.404724821517878  -2.338041431261778   2.507602188919887 ;
%!       -5.554133527986547   1.687202056051516  -4.161958568738229 ];
%! Be = [-0.200000000000000  -1.250000000000000 ;
%!                        0  -0.609718066426859 ;
%!        0.000000000000000   2.221653411194620 ];
%! Ce = [                 0  -0.867926073205492   0.211918650529080 ;
%!                        0                   0   0.900230427447532 ];
%! De = [ 1 0 ;
%!        0 1 ];
%!
%! # State space representations have different structures, therefore
%! # compare the Markov parameters instead
%! M = zeros (2,20);
%! Me = zeros (2,20);
%! for i = 0:9, M(:,2*i+1:2*i+2) = C*A^i*B; Me(:,2*i+1:2*i+2) = Ce*Ae^i*Be; end;
%!
%! assert (M, Me, 1e-6);

