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
## @deftypefn {Function File} {@var{sys} =} minreal (@var{sys})
## @deftypefnx {Function File} {@var{sys} =} minreal (@var{sys}, @var{tol})
## Minimal realization or zero-pole cancellation of @acronym{LTI} models.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.4

function sys = minreal (sys, tol = "def")

  if (nargin > 2)  # nargin == 0 not possible because minreal is inside @lti
    print_usage ();
  endif

  if (! is_real_scalar (tol) && tol != "def")
    error ("minreal: second argument must be a real-valued scalar");
  endif

  sys = __minreal__ (sys, tol);

endfunction


## ss: minreal (SLICOT TB01PD)
%!shared C, D
%!
%! A = ss (-2, 3, 4, 5);
%! B = A / A;
%! C = minreal (B, 1e-15);
%! D = ss (1);
%!
%!assert (C.a, D.a);
%!assert (C.b, D.b);
%!assert (C.c, D.c);
%!assert (C.d, D.d);

%!shared M, Me
%! A = [ 1.0   2.0   0.0
%!       4.0  -1.0   0.0
%!       0.0   0.0   1.0 ];
%!
%! B = [ 1.0
%!       0.0
%!       1.0 ];
%!
%! C = [ 0.0   1.0  -1.0
%!       0.0   0.0   1.0 ];
%!
%! D = zeros (2, 1);
%!
%! [Ar, Br, Cr] = __sl_tb01pd__ (A, B, C, 0.0, true);
%! M = [Ar, Br; Cr, D];
%!
%! Ae = [ 1.0000  -1.4142   1.4142
%!       -2.8284   0.0000   1.0000
%!        2.8284   1.0000   0.0000 ];
%!
%! Be = [-1.0000
%!        0.7071
%!        0.7071 ];
%!
%! Ce = [ 0.0000   0.0000  -1.4142
%!        0.0000   0.7071   0.7071 ];
%!
%! De = zeros (2, 1);
%!
%! Me = [Ae, Be; Ce, De];
%!
%!assert (M, Me, 1e-4);


## dss: minreal (SLICOT TG01JD)
## FIXME: Test fails with larger ldwork in sltg01jd.cc
%!shared Ar, Br, Cr, Dr, Er, Ae, Be, Ce, De, Ee,  num, den, num2, den2, num3, num3a, den3, den3a
%! A = [ -2    -3     0     0     0     0     0     0     0
%!        1     0     0     0     0     0     0     0     0
%!        0     0    -2    -3     0     0     0     0     0
%!        0     0     1     0     0     0     0     0     0
%!        0     0     0     0     1     0     0     0     0
%!        0     0     0     0     0     1     0     0     0
%!        0     0     0     0     0     0     1     0     0
%!        0     0     0     0     0     0     0     1     0
%!        0     0     0     0     0     0     0     0     1 ];
%!
%! E = [  1     0     0     0     0     0     0     0     0
%!        0     1     0     0     0     0     0     0     0
%!        0     0     1     0     0     0     0     0     0
%!        0     0     0     1     0     0     0     0     0
%!        0     0     0     0     0     0     0     0     0
%!        0     0     0     0     1     0     0     0     0
%!        0     0     0     0     0     0     0     0     0
%!        0     0     0     0     0     0     1     0     0
%!        0     0     0     0     0     0     0     1     0 ];
%!
%! B = [  1     0
%!        0     0
%!        0     1
%!        0     0
%!       -1     0
%!        0     0
%!        0    -1
%!        0     0
%!        0     0 ];
%!
%! C = [  1     0     1    -3     0     1     0     2     0
%!        0     1     1     3     0     1     0     0     1 ];
%!
%! D = zeros (2, 2);
%!
%! sys = dss (A, B, C, D, E, "scaled", true);
%! sysmin = minreal (sys, 0.0);
%! [Ar, Br, Cr, Dr, Er] = dssdata (sysmin);
%! [num,den]=tfdata(sys);
%! sysmin = minreal (sys, 1e-6);
%! [num2,den2]=tfdata(sysmin);
%! sys3 = dss (Ar, Br, Cr, Dr, Er, "scaled", true);
%! [num3,den3]=tfdata(sys3);
%! sysmin3 = minreal (sys3, 1e-6);
%! [num3a,den3a]=tfdata(sysmin3);
%!
%!assert (num, num2 , 1e-4);
%!assert (den, den2 , 1e-4);
%!assert (num3, num3a , 1e-4);
%!assert (den3, den3a , 1e-4);
%!assert (num, num3a, 1e-4);
%!assert (den, den3a, 1e-4);


## tf: minreal
%!shared a, b, c, d
%! s = tf ("s");
%! G1 = (s+1)*s*5/(s+1)/(s^2+s+1);
%! G2 = tf ([1, 1, 1], [2, 2, 2]);
%! G1min = minreal (G1);
%! G2min = minreal (G2);
%! a = G1min.num{1, 1};
%! b = G1min.den{1, 1};
%! c = G2min.num{1, 1};
%! d = G2min.den{1, 1};
%!assert (a, [5, 0], 1e-4);
%!assert (b, [1, 1, 1], 1e-4);
%!assert (c, 0.5, 1e-4);
%!assert (d, 1, 1e-4);
