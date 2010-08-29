## Copyright (C) 2009   Luca Favatella
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
## @deftypefn {Function File} test ltimodels
## @deftypefnx {Function File} ltimodels
## @deftypefnx {Function File} ltimodels (@var{systype})
## Test suite and help for LTI models.
## @end deftypefn

## Author: Luca Favatella <slackydeb@gmail.com>
## Created: November 2009
## Version: 0.1

function ltimodels (systype = "general")

  %if (nargin > 1)
    print_usage ();
  %endif

  if (! ischar (systype))
    error ("ltimodels: argument must be a string");
  endif

  systype = lower (systype);

  switch (systype)
    case "ss"
      str = {"State Space (SS) Models"\
             "-----------------------"\
             ""};

    case "tf"
      str = {"Transfer Function (TF) Models"\
             "-----------------------------"\
             ""};

    otherwise  # general
      str = {"Linear Time Invariant (LTI) Models"\
             "----------------------------------"\
             ""};

  endswitch

  disp ("");
  disp (char (str));

endfunction


## isct, isdt
%!shared ltisys
%! ltisys = tf (12);
%!assert (ltisys.ts, -1);
%!assert (isct (ltisys));
%!assert (isdt (ltisys));

%!shared ltisys
%! ltisys = ss (17);
%!assert (ltisys.ts, -1);
%!assert (isct (ltisys));
%!assert (isdt (ltisys));

%!shared ltisys
%! ltisys = tf (1, [1 1]);
%!assert (ltisys.ts, 0);
%!assert (isct (ltisys));
%!assert (! isdt (ltisys));

%!shared ltisys, ts
%! ts = 0.1;
%! ltisys = ss (-1, 1, 1, 0, ts);
%!assert (ltisys.ts, ts);
%!assert (! isct (ltisys));
%!assert (isdt (ltisys));

## mtimes
%!shared sysmat, sysmat_exp
%! sys1 = ss ([0, 1; -3, -2], [0; 1], [-5, 1], [2]);
%! sys2 = ss ([-10], [1], [-40], [5]);
%! sys3 = sys2 * sys1;
%! [A, B, C, D] = ssdata (sys3);
%! sysmat = [A, B; C, D];
%! A_exp = [ -10   -5    1
%!             0    0    1
%!             0   -3   -2 ];
%! B_exp = [   2
%!             0
%!             1 ];
%! C_exp = [ -40  -25    5 ];
%! D_exp = [  10 ];
%! sysmat_exp = [A_exp, B_exp; C_exp, D_exp];
%!assert (sysmat, sysmat_exp)

## norm ct
%!shared H2, Hinf
%! sys = ss (-1, 1, 1, 0);
%! H2 = norm (sys, 2);
%! Hinf = norm (sys, inf);
%!assert (H2, 0.7071, 1.5e-5);
%!assert (Hinf, 1, 5e-4);

## norm dt
%!shared H2, Hinf
%! a = [ 2.417   -1.002    0.5488
%!           2        0         0
%!           0      0.5         0 ];
%! b = [     1
%!           0
%!           0 ];
%! c = [-0.424    0.436   -0.4552 ];
%! d = [     1 ];
%! sys = ss (a, b, c, d, 0.1);
%! H2 = norm (sys, 2);
%! Hinf = norm (sys, inf);
%!assert (H2, 1.2527, 1.5e-5);
%!assert (Hinf, 2.7, 0.1);

## transmission zeros of state-space models

## Results from the "Dark Side" 7.5 and 7.8
##
##  -13.2759
##   12.5774
##  -0.0155

## Results from Scilab 5.2.0b1 (trzeros)
##
##  - 13.275931  
##    12.577369  
##  - 0.0155265

%!shared z, z_exp
%! A = [   -0.7   -0.0458     -12.2        0
%!            0    -0.014   -0.2904   -0.562
%!            1   -0.0057      -1.4        0
%!            1         0         0        0 ];
%!
%! B = [  -19.1      -3.1
%!      -0.0119   -0.0096
%!        -0.14     -0.72
%!            0         0 ];
%!
%! C = [      0         0        -1        1
%!            0         0     0.733        0 ];
%!
%! D = [      0         0
%!       0.0768    0.1134 ];
%!
%! sys = ss (A, B, C, D);
%! z = sort (zero (sys));
%!
%! z_exp = sort ([-13.2759; 12.5774; -0.0155]);
%!
%!assert (z, z_exp, 1e-4);

## inverse of state-space models
## test from SLICOT AB07ND
## note the negative signs in Me for compatibility reasons
%!shared M, Me
%! A = [ 1.0   2.0   0.0
%!       4.0  -1.0   0.0
%!       0.0   0.0   1.0 ];
%!
%! B = [ 1.0   0.0
%!       0.0   1.0
%!       1.0   0.0 ];
%!
%! C = [ 0.0   1.0  -1.0
%!       0.0   0.0   1.0 ];
%!
%! D = [ 4.0   0.0
%!       0.0   1.0 ];
%!
%! sys = ss (A, B, C, D);
%! sysinv = inv (sys);
%! [Ai, Bi, Ci, Di] = ssdata (sysinv);
%! M = [Ai, Bi; Ci, Di];
%!
%! Ae = [ 1.0000   1.7500   0.2500
%!        4.0000  -1.0000  -1.0000
%!        0.0000  -0.2500   1.2500 ];
%!
%! Be = [-0.2500   0.0000
%!        0.0000  -1.0000
%!       -0.2500   0.0000 ];
%!
%! Ce = [ 0.0000   0.2500  -0.2500
%!        0.0000   0.0000   1.0000 ];
%!
%! De = [ 0.2500   0.0000
%!        0.0000   1.0000 ];
%!
%! Me = [Ae, -Be; -Ce, De];
%!
%!assert (M, Me, 1e-4);

## ss: minreal
%!shared C, D
%!
%! A = ss (-2, 3, 4, 5);
%! B = A / A;
%! C = minreal (B);
%! D = ss (1);
%!
%!assert (C.a, D.a);
%!assert (C.b, D.b);
%!assert (C.c, D.c);
%!assert (C.d, D.d);

