## Copyright (C) 2011   Lukas F. Reichlin
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
## @deftypefn{Function File} {[@var{K}, @var{N}, @var{gamma}] =} ncfsyn (@var{G}, @var{W1}, @var{W2}, @var{factor})
## Compute positive feedback controller using the McFarlane/Glover Loop Shaping Design Procedure.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: July 2011
## Version: 0.1

% function [K, N, gamma] = ncfsyn (G, W1 = [], W2 = [], factor = 1.0)
function K = ncfsyn (G, W1 = [], W2 = [], factor = 1.0)

  if (nargin == 0 || nargin > 4)
    print_usage ();
  endif
  
  if (! isa (G, "lti"))
    error ("ncfsyn: first argument must be an LTI system");
  endif
  
  if (! is_real_scalar (factor) || factor < 1.0)
    error ("ncfsyn: fourth argument invalid");
  endif

  [p, m] = size (G);

  W1 = __adjust_weighting__ (W1, m);
  W2 = __adjust_weighting__ (W2, p);
  
  Gs = W2 * G * W1;            # shaped plant

  [a, b, c, d, tsam] = ssdata (Gs);

  ## synthesis
  if (isct (Gs))               # continuous-time
    [ak, bk, ck, dk, rcond] = slsb10id (a, b, c, d, factor);
  elseif (any (d(:)))          # discrete-time, d != 0
    [ak, bk, ck, dk, rcond] = slsb10zd (a, b, c, d, factor);
  else                         # discrete-time, d == 0
    [ak, bk, ck, dk, rcond] = slsb10kd (a, b, c, d, factor);
  endif

  ## controller
  K = ss (ak, bk, ck, dk, tsam);

endfunction


function W = __adjust_weighting__ (W, s)

  if (isempty (W))
    W = ss (eye (s));
  else
    W = ss (W);
    [p, m] = size (W);
    if (m == s && p == s)      # model is of correct size
      return;
    elseif (m == 1 && p == 1)  # model is SISO
      tmp = W;
      for k = 2 : s
        W = append (W, tmp);   # stack SISO model s times
      endfor
    else                       # model is invalid
      error ("ncfsyn: %s must have 1 or %d inputs", inputname (1), s);
    endif
  endif

endfunction


## continuous-time case, direct access to sb10id
%!shared AK, BK, CK, DK, RCOND, AKe, BKe, CKe, DKe, RCONDe
%! A = [ -1.0  0.0  4.0  5.0 -3.0 -2.0
%!       -2.0  4.0 -7.0 -2.0  0.0  3.0
%!       -6.0  9.0 -5.0  0.0  2.0 -1.0
%!       -8.0  4.0  7.0 -1.0 -3.0  0.0
%!        2.0  5.0  8.0 -9.0  1.0 -4.0
%!        3.0 -5.0  8.0  0.0  2.0 -6.0 ];
%!   
%! B = [ -3.0 -4.0
%!        2.0  0.0
%!       -5.0 -7.0
%!        4.0 -6.0
%!       -3.0  9.0
%!        1.0 -2.0 ];
%!   
%! C = [  1.0 -1.0  2.0 -4.0  0.0 -3.0
%!       -3.0  0.0  5.0 -1.0  1.0  1.0
%!       -7.0  5.0  0.0 -8.0  2.0 -2.0 ];
%!  
%! D = [  1.0 -2.0
%!        0.0  4.0
%!        5.0 -3.0 ];
%!    
%! FACTOR = 1.0;
%!
%! [AK, BK, CK, DK, RCOND] = slsb10id (A, B, C, D, FACTOR);
%!
%! AKe = [ -39.0671    9.9293   22.2322  -27.4113   43.8655
%!          -6.6117    3.0006   11.0878  -11.4130   15.4269
%!          33.6805   -6.6934  -23.9953   14.1438  -33.4358
%!         -32.3191    9.7316   25.4033  -24.0473   42.0517
%!         -44.1655   18.7767   34.8873  -42.4369   50.8437 ];
%!
%! BKe = [ -10.2905  -16.5382  -10.9782
%!          -4.3598   -8.7525   -5.1447
%!           6.5962    1.8975    6.2316
%!          -9.8770  -14.7041  -11.8778
%!          -9.6726  -22.7309  -18.2692 ];
%!
%! CKe = [  -0.6647   -0.0599   -1.0376    0.5619    1.7297
%!          -8.4202    3.9573    7.3094   -7.6283   10.6768 ];
%!
%! DKe = [  0.8466    0.4979   -0.6993
%!         -1.2226   -4.8689   -4.5056 ];
%!
%! RCONDe = [ 0.13861D-01  0.90541D-02 ].';
%!
%!assert (AK, AKe, 1e-4);
%!assert (BK, BKe, 1e-4);
%!assert (CK, CKe, 1e-4);
%!assert (DK, DKe, 1e-4);
%!assert (RCOND, RCONDe, 1e-4);


## continuous-time case
%!shared AK, BK, CK, DK, RCOND, AKe, BKe, CKe, DKe, RCONDe
%! A = [ -1.0  0.0  4.0  5.0 -3.0 -2.0
%!       -2.0  4.0 -7.0 -2.0  0.0  3.0
%!       -6.0  9.0 -5.0  0.0  2.0 -1.0
%!       -8.0  4.0  7.0 -1.0 -3.0  0.0
%!        2.0  5.0  8.0 -9.0  1.0 -4.0
%!        3.0 -5.0  8.0  0.0  2.0 -6.0 ];
%!   
%! B = [ -3.0 -4.0
%!        2.0  0.0
%!       -5.0 -7.0
%!        4.0 -6.0
%!       -3.0  9.0
%!        1.0 -2.0 ];
%!   
%! C = [  1.0 -1.0  2.0 -4.0  0.0 -3.0
%!       -3.0  0.0  5.0 -1.0  1.0  1.0
%!       -7.0  5.0  0.0 -8.0  2.0 -2.0 ];
%!  
%! D = [  1.0 -2.0
%!        0.0  4.0
%!        5.0 -3.0 ];
%!    
%! FACTOR = 1.0;
%!
%! G = ss (A, B, C, D);
%! K = ncfsyn (G, [], [], FACTOR);
%! [AK, BK, CK, DK] = ssdata (K);
%!
%! AKe = [ -39.0671    9.9293   22.2322  -27.4113   43.8655
%!          -6.6117    3.0006   11.0878  -11.4130   15.4269
%!          33.6805   -6.6934  -23.9953   14.1438  -33.4358
%!         -32.3191    9.7316   25.4033  -24.0473   42.0517
%!         -44.1655   18.7767   34.8873  -42.4369   50.8437 ];
%!
%! BKe = [ -10.2905  -16.5382  -10.9782
%!          -4.3598   -8.7525   -5.1447
%!           6.5962    1.8975    6.2316
%!          -9.8770  -14.7041  -11.8778
%!          -9.6726  -22.7309  -18.2692 ];
%!
%! CKe = [  -0.6647   -0.0599   -1.0376    0.5619    1.7297
%!          -8.4202    3.9573    7.3094   -7.6283   10.6768 ];
%!
%! DKe = [  0.8466    0.4979   -0.6993
%!         -1.2226   -4.8689   -4.5056 ];
%!
%! RCONDe = [ 0.13861D-01  0.90541D-02 ];
%!
%!assert (AK, AKe, 1e-4);
%!assert (BK, BKe, 1e-4);
%!assert (CK, CKe, 1e-4);
%!assert (DK, DKe, 1e-4);
