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
## @deftypefn {Function File} {@var{sys} =} c2d (@var{sys}, @var{tsam})
## @deftypefnx {Function File} {@var{sys} =} c2d (@var{sys}, @var{tsam}, @var{method})
## @deftypefnx {Function File} {@var{sys} =} c2d (@var{sys}, @var{tsam}, @var{'prewarp'}, @var{w0})
## Convert the continuous @acronym{LTI} model into its discrete-time equivalent.
##
## @strong{Inputs}
## @table @var
## @item sys
## Continuous-time @acronym{LTI} model.
## @item tsam
## Sampling time in seconds.
## @item method
## Optional conversion method.  If not specified, default method @var{"zoh"}
## is taken.
## @table @var
## @item 'impulse'
## Impulse Invarient transformation.
## @item 'zoh'
## Zero-order hold or matrix exponential.
## @item 'tustin', 'bilin'
## Bilinear transformation or Tustin approximation.
## @item 'prewarp'
## Bilinear transformation with pre-warping at frequency @var{w0}.
## @item 'matched'
## Matched pole/zero method.
## @end table
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## Discrete-time @acronym{LTI} model.
## @end table
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.3

function sys = c2d (sys, tsam, method = "std", w0 = 0)

  if (nargin < 2 || nargin > 4)
    print_usage ();
  endif

  if (! isa (sys, "lti"))
    error ("c2d: first argument is not an LTI model");
  endif

  if (isdt (sys))
    error ("c2d: system is already discrete-time");
  endif

  if (! issample (tsam))
    error ("c2d: second argument is not a valid sample time");
  endif

  if (! ischar (method))
    error ("c2d: third argument is not a string");
  endif

  if (! issample (w0, 0))
    error ("c2d: fourth argument is not a valid pre-warping frequency");
  endif

  sys = __c2d__ (sys, tsam, lower (method), w0);
  sys.tsam = tsam;

endfunction


## bilinear transformation
## using oct-file directly
%!shared Mo, Me
%! A = [  1.0  0.5
%!        0.5  1.0 ].';
%!
%! B = [  0.0 -1.0
%!        1.0  0.0 ].';
%!
%! C = [ -1.0  0.0
%!        0.0  1.0 ].';
%!
%! D = [  1.0  0.0
%!        0.0 -1.0 ].';
%!
%! [Ao, Bo, Co, Do] = __sl_ab04md__ (A, B, C, D, 1.0, 1.0, false);
%!
%! Ae = [ -1.0000  -4.0000
%!        -4.0000  -1.0000 ];
%!
%! Be = [  2.8284   0.0000
%!         0.0000  -2.8284 ];
%!
%! Ce = [  0.0000   2.8284
%!        -2.8284   0.0000 ];
%!
%! De = [ -1.0000   0.0000
%!         0.0000  -3.0000 ];
%!
%! Mo = [Ao, Bo; Co, Do];
%! Me = [Ae, Be; Ce, De];
%!
%!assert (Mo, Me, 1e-4);


## bilinear transformation
## user function
%!shared Mo, Me
%! A = [  1.0  0.5
%!        0.5  1.0 ].';
%!
%! B = [  0.0 -1.0
%!        1.0  0.0 ].';
%!
%! C = [ -1.0  0.0
%!        0.0  1.0 ].';
%!
%! D = [  1.0  0.0
%!        0.0 -1.0 ].';
%!
%! [Ao, Bo, Co, Do] = ssdata (c2d (ss (A, B, C, D), 2, "tustin"));
%!
%! Ae = [ -1.0000  -4.0000
%!        -4.0000  -1.0000 ];
%!
%! Be = [  2.8284   0.0000
%!         0.0000  -2.8284 ];
%!
%! Ce = [  0.0000   2.8284
%!        -2.8284   0.0000 ];
%!
%! De = [ -1.0000   0.0000
%!         0.0000  -3.0000 ];
%!
%! Mo = [Ao, Bo; Co, Do];
%! Me = [Ae, Be; Ce, De];
%!
%!assert (Mo, Me, 1e-4);


## bilinear transformation
## both directions
%!shared Mo, Me
%! A = [  1.0  0.5
%!        0.5  1.0 ];
%!
%! B = [  0.0 -1.0
%!        1.0  0.0 ];
%!
%! C = [ -1.0  0.0
%!        0.0  1.0 ];
%!
%! D = [  1.0  0.0
%!        0.0 -1.0 ];
%!
%! [Ao, Bo, Co, Do] = ssdata (d2c (c2d (ss (A, B, C, D), 2, "tustin"), "tustin"));
%!
%! Mo = [Ao, Bo; Co, Do];
%! Me = [A, B; C, D];
%!
%!assert (Mo, Me, 1e-4);


## zero-order hold
## both directions
%!shared Mo, Me
%! A = [  1.0  0.5
%!        0.5  1.0 ];
%!
%! B = [  0.0 -1.0
%!        1.0  0.0 ];
%!
%! C = [ -1.0  0.0
%!        0.0  1.0 ];
%!
%! D = [  1.0  0.0
%!        0.0 -1.0 ];
%!
%! [Ao, Bo, Co, Do] = ssdata (d2c (c2d (ss (A, B, C, D), 2, "zoh"), "zoh"));
%!
%! Mo = [Ao, Bo; Co, Do];
%! Me = [A, B; C, D];
%!
%!assert (Mo, Me, 1e-4);


## bilinear transformation with pre-warping
## both directions
%!shared Mo, Me
%! A = [  1.0  0.5
%!        0.5  1.0 ];
%!
%! B = [  0.0 -1.0
%!        1.0  0.0 ];
%!
%! C = [ -1.0  0.0
%!        0.0  1.0 ];
%!
%! D = [  1.0  0.0
%!        0.0 -1.0 ];
%!
%! [Ao, Bo, Co, Do] = ssdata (d2c (c2d (ss (A, B, C, D), 2, "prewarp", 1000), "prewarp", 1000));
%!
%! Mo = [Ao, Bo; Co, Do];
%! Me = [A, B; C, D];
%!
%!assert (Mo, Me, 1e-4);

## matrix exponential
%!shared Aex, Aexint, Aex_exp, Aexint_exp
%! A = [  5.0   4.0   3.0   2.0   1.0
%!        1.0   6.0   0.0   4.0   3.0
%!        2.0   0.0   7.0   6.0   5.0
%!        1.0   3.0   1.0   8.0   7.0
%!        2.0   5.0   7.0   1.0   9.0 ];
%!
%! Aex_exp = [ 1.8391   0.9476   0.7920   0.8216   0.7811
%!             0.3359   2.2262   0.4013   1.0078   1.0957
%!             0.6335   0.6776   2.6933   1.6155   1.8502
%!             0.4804   1.1561   0.9110   2.7461   2.0854
%!             0.7105   1.4244   1.8835   1.0966   3.4134 ];
%!
%! Aexint_exp = [ 0.1347   0.0352   0.0284   0.0272   0.0231
%!                0.0114   0.1477   0.0104   0.0369   0.0368
%!                0.0218   0.0178   0.1624   0.0580   0.0619
%!                0.0152   0.0385   0.0267   0.1660   0.0732
%!                0.0240   0.0503   0.0679   0.0317   0.1863 ];
%!
%! [Aex, Aexint] = __sl_mb05nd__ (A, 0.1, 0.0001);
%!
%!assert (Aex, Aex_exp, 1e-4);
%!assert (Aexint, Aexint_exp, 1e-4);


