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
## @deftypefn{Function File} {[@var{K}, @var{N}, @var{gamma}, @var{info}] =} hinfsyn (@var{P}, @var{nmeas}, @var{ncon})
## @deftypefnx{Function File} {[@var{K}, @var{N}, @var{gamma}, @var{info}] =} hinfsyn (@var{P}, @var{nmeas}, @var{ncon}, @dots{})
## @deftypefnx{Function File} {[@var{K}, @var{N}, @var{gamma}, @var{info}] =} hinfsyn (@var{P}, @var{nmeas}, @var{ncon}, @var{opt}, @dots{})
## @deftypefnx{Function File} {[@var{K}, @var{N}, @var{gamma}, @var{info}] =} hinfsyn (@var{P}, @dots{})
## @deftypefnx{Function File} {[@var{K}, @var{N}, @var{gamma}, @var{info}] =} hinfsyn (@var{P}, @var{opt}, @dots{})
## H-infinity control synthesis for @acronym{LTI} plant.
##
## @strong{Inputs}
## @table @var
## @item P
## Generalized plant.  Must be a proper/realizable @acronym{LTI} model.
## If @var{P} is constructed with @command{mktito} or @command{augw},
## arguments @var{nmeas} and @var{ncon} can be omitted.
## @item nmeas
## Number of measured outputs v.  The last @var{nmeas} outputs of @var{P} are connected to the
## inputs of controller @var{K}.  The remaining outputs z (indices 1 to p-nmeas) are used
## to calculate the H-infinity norm.
## @item ncon
## Number of controlled inputs u.  The last @var{ncon} inputs of @var{P} are connected to the
## outputs of controller @var{K}.  The remaining inputs w (indices 1 to m-ncon) are excited
## by a harmonic test signal.
## @item @dots{}
## Optional pairs of keys and values.  @code{'key1', value1, 'key2', value2}.
## @item opt
## Optional struct with keys as field names.
## Struct @var{opt} can be created directly or
## by function @command{options}.  @code{opt.key1 = value1, opt.key2 = value2}.
## @end table
##
## @strong{Outputs}
## @table @var
## @item K
## State-space model of the H-infinity (sub-)optimal controller.
## @item N
## State-space model of the lower LFT of @var{P} and @var{K}.
## @item info
## Structure containing additional information.
## @item info.gamma
## L-infinity norm of @var{N}.
## @item info.rcond
## Vector @var{rcond} contains estimates of the reciprocal condition
## numbers of the matrices which are to be inverted and
## estimates of the reciprocal condition numbers of the
## Riccati equations which have to be solved during the
## computation of the controller @var{K}.  For details,
## see the description of the corresponding SLICOT routine.
## @end table
##
## @strong{Option Keys and Values}
## @table @var
## @item 'method'
## String specifying the desired kind of controller:
## @table @var
## @item 'optimal', 'opt', 'o'
## Compute optimal controller using gamma iteration.
## Default selection for compatibility reasons.
## @item 'suboptimal', 'sub', 's'
## Compute (sub-)optimal controller.  For stability reasons,
## suboptimal controllers are to be preferred over optimal ones.
## @end table
## @item 'gmax'
## The maximum value of the H-infinity norm of @var{N}.
## It is assumed that @var{gmax} is sufficiently large
## so that the controller is admissible.  Default value is 1e15.
## @item 'gmin'
## Initial lower bound for gamma iteration.  Default value is 0.
## @var{gmin} is only meaningful for optimal discrete-time controllers.
## @item 'tolgam'
## Tolerance used for controlling the accuracy of @var{gamma}
## and its distance to the estimated minimal possible
## value of @var{gamma}.  Default value is 0.01.
## If @var{tolgam} = 0, then a default value equal to @code{sqrt(eps)}
## is used, where @var{eps} is the relative machine precision.
## For suboptimal controllers, @var{tolgam} is ignored.
## @item 'actol'
## Upper bound for the poles of the closed-loop system @var{N}
## used for determining if it is stable.
## @var{actol} >= 0 for stable systems.
## For suboptimal controllers, @var{actol} is ignored.
## @end table
##
## @strong{Block Diagram}
## @example
## @group
##
## gamma = min||N(K)||             N = lft (P, K)
##          K         inf
##
##                +--------+  
##        w ----->|        |-----> z
##                |  P(s)  |
##        u +---->|        |-----+ v
##          |     +--------+     |
##          |                    |
##          |     +--------+     |
##          +-----|  K(s)  |<----+
##                +--------+
##
##                +--------+      
##        w ----->|  N(s)  |-----> z
##                +--------+
## @end group
## @end example
##
## @strong{Algorithm}@*
## Uses SLICOT SB10FD, SB10DD and SB10AD by courtesy of
## @uref{http://www.slicot.org, NICONET e.V.}
##
## @seealso{augw, mixsyn}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: December 2009
## Version: 0.3

function [K, varargout] = hinfsyn (P, varargin)

  ## check input arguments
  if (nargin == 0)
    print_usage ();
  endif

  if (! isa (P, "lti"))
    error ("hinfsyn: first argument must be an LTI system");
  endif
  
  if (nargin == 1 || (nargin > 1 && ! is_real_scalar (varargin{1})))    # hinfsyn (P, ...)
    [nmeas, ncon] = __tito_dim__ (P, "hinfsyn");
  elseif (nargin >= 3)                          # hinfsyn (P, nmeas, ncon, ...)
    nmeas = varargin{1};
    ncon = varargin{2};
    varargin = varargin(3:end);
  else
    print_usage ();
  endif
  
  if (! is_real_scalar (nmeas))
    error ("hinfsyn: second argument 'nmeas' invalid");
  endif
  
  if (! is_real_scalar (ncon))
    error ("hinfsyn: third argument 'ncon' invalid");
  endif
  
  if (numel (varargin) > 0 && isstruct (varargin{1}))   # hinfsyn (P, nmeas, ncon, opt, ...), hinfsyn (P, opt, ...)
    varargin = horzcat (__opt2cell__ (varargin{1}), varargin(2:end));
  endif

  nkv = numel (varargin);   # number of keys and values
  
  if (rem (nkv, 2))
    error ("hinfsyn: keys and values must come in pairs");
  endif
  
  ## default arguments
  gmax = 1e15;
  gmin = 0;
  tolgam = 0.01;
  actol = eps;              # tolerance for stability margin
  method = "opt";
  
  ## handle keys and values
  for k = 1 : 2 : nkv
    key = lower (varargin{k});
    val = varargin{k+1};
    switch (key)
      case "gmax"
        if (! is_real_scalar (val) || val < 0)
          error ("hinfsyn: 'gmax' must be a real-valued, non-negative scalar");
        endif
        gmax = val;
      case "gmin"
        if (! is_real_scalar (val) || val < 0)
          error ("hinfsyn: 'gmin' must be a real-valued, non-negative scalar");
        endif
        gmin = val;
      case "tolgam"    
        if (! is_real_scalar (val) || val < 0)
          error ("hinfsyn: 'tolgam' must be a real-valued, non-negative scalar");
        endif
        tolgam = val;
      case "actol"  
        if (! is_real_scalar (val) || val < 0)
          error ("hinfsyn: 'actol' must be a real-valued, non-negative scalar");
        endif
        actol = val;
      case "method"
        ## NOTE: I called this "method" because of the dark side,
        ##       maybe something like "type" would make more sense ...               
        if (strncmpi (val, "s", 1))
          method = "sub";   # sub-optimal
        elseif (strncmpi (val, "o", 1) || strncmpi (val, "ric", 1))
          method = "opt";   # optimal
        else
          error ("hinfsyn: invalid method '%s'", val);
        endif
      otherwise
        warning ("hinfsyn: invalid property name '%s' ignored", key);
    endswitch
  endfor

  [a, b, c, d, tsam] = ssdata (P);
  
  ## check assumption A1
  m = columns (b);
  p = rows (c);
  
  m1 = m - ncon;
  p1 = p - nmeas;
  
  if (! isstabilizable (P(:, m1+1:m)))
    error ("hinfsyn: (A, B2) must be stabilizable");
  endif
  
  if (! isdetectable (P(p1+1:p, :)))
    error ("hinfsyn: (C2, A) must be detectable");
  endif

  ## H-infinity synthesis
  switch (method)
    case "sub"              # sub-optimal controller
      if (isct (P))         # continuous-time plant
        [ak, bk, ck, dk, rcond] = __sl_sb10fd__ (a, b, c, d, ncon, nmeas, gmax);
      else                  # discrete-time plant
        [ak, bk, ck, dk, rcond] = __sl_sb10dd__ (a, b, c, d, ncon, nmeas, gmax);
      endif

    case "opt"              # optimal controller
      if (isct (P))         # continuous-time plant
        [ak, bk, ck, dk, ~, ~, ~, ~, ~, rcond] = __sl_sb10ad__ (a, b, c, d, ncon, nmeas, gmax, tolgam, -actol);
      else                  # discrete-time plant
        ## NOTE: check whether it is an alternative to compute the bilinear transformation
        ##       of P, use __sl_sb10ad__ for a continuous-time controller and then
        ##       discretize the controller.
        ## estimate gamma
        Pt = d2c (P, "tustin");
        [at, bt, ct, dt] = ssdata (Pt);
        [~, ~, ~, ~, ~, ~, ~, ~, gamma] = __sl_sb10ad__ (at, bt, ct, dt, ncon, nmeas, gmax, tolgam, -actol);
        ## gamma iteration - bisection method using __sl_sb10dd__
        gmax = 1.2*gamma;
        while (gmax > eps && (gmax - gmin)/gmax > tolgam)
          gmid = (gmax + gmin)/2;
          try
            [ak, bk, ck, dk, rcond] = __sl_sb10dd__ (a, b, c, d, ncon, nmeas, gmid);
            ## check for stability
            K = ss (ak, bk, ck, dk, tsam);
            N = lft (P, K);
            if (isstable (N, actol))
              gmax = norm (N, inf);
            else
              gmin = gmid;
            endif
          catch             # cannot find solution
            gmin = gmid;
          end_try_catch 
        endwhile
      endif
    
    otherwise
      error ("hinfsyn: this should never happen");
  endswitch
  
  ## controller
  K = ss (ak, bk, ck, dk, tsam);
  
  if (nargout > 1)
    N = lft (P, K);
    varargout{1} = N;
    if (nargout > 2)
      gamma = norm (N, inf);
      varargout{2} = gamma;
      if (nargout > 3)
        varargout{3} = struct ("gamma", gamma, "rcond", rcond);
      endif
    endif
  endif

endfunction


## sub-optimal controller, continuous-time case
%!shared M, M_exp
%! A = [-1.0  0.0  4.0  5.0 -3.0 -2.0
%!      -2.0  4.0 -7.0 -2.0  0.0  3.0
%!      -6.0  9.0 -5.0  0.0  2.0 -1.0
%!      -8.0  4.0  7.0 -1.0 -3.0  0.0
%!       2.0  5.0  8.0 -9.0  1.0 -4.0
%!       3.0 -5.0  8.0  0.0  2.0 -6.0];
%!
%! B = [-3.0 -4.0 -2.0  1.0  0.0
%!       2.0  0.0  1.0 -5.0  2.0
%!      -5.0 -7.0  0.0  7.0 -2.0
%!       4.0 -6.0  1.0  1.0 -2.0
%!      -3.0  9.0 -8.0  0.0  5.0
%!       1.0 -2.0  3.0 -6.0 -2.0];
%!
%! C = [ 1.0 -1.0  2.0 -4.0  0.0 -3.0
%!      -3.0  0.0  5.0 -1.0  1.0  1.0
%!      -7.0  5.0  0.0 -8.0  2.0 -2.0
%!       9.0 -3.0  4.0  0.0  3.0  7.0
%!       0.0  1.0 -2.0  1.0 -6.0 -2.0];
%!
%! D = [ 1.0 -2.0 -3.0  0.0  0.0
%!       0.0  4.0  0.0  1.0  0.0
%!       5.0 -3.0 -4.0  0.0  1.0
%!       0.0  1.0  0.0  1.0 -3.0
%!       0.0  0.0  1.0  7.0  1.0];
%!
%! P = ss (A, B, C, D);
%! K = hinfsyn (P, 2, 2, "method", "sub", "gmax", 15);
%! M = [K.A, K.B; K.C, K.D];
%!
%! KA = [ -2.8043  14.7367   4.6658   8.1596   0.0848   2.5290
%!         4.6609   3.2756  -3.5754  -2.8941   0.2393   8.2920
%!       -15.3127  23.5592  -7.1229   2.7599   5.9775  -2.0285
%!       -22.0691  16.4758  12.5523 -16.3602   4.4300  -3.3168
%!        30.6789  -3.9026  -1.3868  26.2357  -8.8267  10.4860
%!        -5.7429   0.0577  10.8216 -11.2275   1.5074 -10.7244];
%!
%! KB = [ -0.1581  -0.0793
%!        -0.9237  -0.5718
%!         0.7984   0.6627
%!         0.1145   0.1496
%!        -0.6743  -0.2376
%!         0.0196  -0.7598];
%!
%! KC = [ -0.2480  -0.1713  -0.0880   0.1534   0.5016  -0.0730
%!         2.8810  -0.3658   1.3007   0.3945   1.2244   2.5690];
%!
%! KD = [  0.0554   0.1334
%!        -0.3195   0.0333];
%!
%! M_exp = [KA, KB; KC, KD];
%!
%!assert (M, M_exp, 1e-4);


## sub-optimal controller, discrete-time case
%!shared M, M_exp
%! A = [-0.7  0.0  0.3  0.0 -0.5 -0.1
%!      -0.6  0.2 -0.4 -0.3  0.0  0.0
%!      -0.5  0.7 -0.1  0.0  0.0 -0.8
%!      -0.7  0.0  0.0 -0.5 -1.0  0.0
%!       0.0  0.3  0.6 -0.9  0.1 -0.4
%!       0.5 -0.8  0.0  0.0  0.2 -0.9];
%!
%! B = [-1.0 -2.0 -2.0  1.0  0.0
%!       1.0  0.0  1.0 -2.0  1.0
%!      -3.0 -4.0  0.0  2.0 -2.0
%!       1.0 -2.0  1.0  0.0 -1.0
%!       0.0  1.0 -2.0  0.0  3.0
%!       1.0  0.0  3.0 -1.0 -2.0];
%!
%! C = [ 1.0 -1.0  2.0 -2.0  0.0 -3.0
%!      -3.0  0.0  1.0 -1.0  1.0  0.0
%!       0.0  2.0  0.0 -4.0  0.0 -2.0
%!       1.0 -3.0  0.0  0.0  3.0  1.0
%!       0.0  1.0 -2.0  1.0  0.0 -2.0];
%!
%! D = [ 1.0 -1.0 -2.0  0.0  0.0
%!       0.0  1.0  0.0  1.0  0.0
%!       2.0 -1.0 -3.0  0.0  1.0
%!       0.0  1.0  0.0  1.0 -1.0
%!       0.0  0.0  1.0  2.0  1.0];
%!
%! P = ss (A, B, C, D, 1);  # value of sampling time doesn't matter
%! K = hinfsyn (P, 2, 2, "method", "sub", "gmax", 111.294);
%! M = [K.A, K.B; K.C, K.D];
%!
%! KA = [-18.0030  52.0376  26.0831  -0.4271 -40.9022  18.0857
%!        18.8203 -57.6244 -29.0938   0.5870  45.3309 -19.8644
%!       -26.5994  77.9693  39.0368  -1.4020 -60.1129  26.6910
%!       -21.4163  62.1719  30.7507  -0.9201 -48.6221  21.8351
%!        -0.8911   4.2787   2.3286  -0.2424  -3.0376   1.2169
%!        -5.3286  16.1955   8.4824  -0.2489 -12.2348   5.1590];
%!
%! KB = [ 16.9788  14.1648
%!       -18.9215 -15.6726
%!        25.2046  21.2848
%!        20.1122  16.8322
%!         1.4104   1.2040
%!         5.3181   4.5149];
%!
%! KC = [ -9.1941  27.5165  13.7364  -0.3639 -21.5983   9.6025
%!         3.6490 -10.6194  -5.2772   0.2432   8.1108  -3.6293];
%!
%! KD = [  9.0317   7.5348
%!        -3.4006  -2.8219];
%!
%! M_exp = [KA, KB; KC, KD];
%!
%!assert (M, M_exp, 1e-4);


## optimal controller, discrete-time case??? -- test for bisection method
%!shared M, M_exp, GAM_exp, GAM
%! A = [-0.7  0.0  0.3  0.0 -0.5 -0.1
%!      -0.6  0.2 -0.4 -0.3  0.0  0.0
%!      -0.5  0.7 -0.1  0.0  0.0 -0.8
%!      -0.7  0.0  0.0 -0.5 -1.0  0.0
%!       0.0  0.3  0.6 -0.9  0.1 -0.4
%!       0.5 -0.8  0.0  0.0  0.2 -0.9];
%!
%! B = [-1.0 -2.0 -2.0  1.0  0.0
%!       1.0  0.0  1.0 -2.0  1.0
%!      -3.0 -4.0  0.0  2.0 -2.0
%!       1.0 -2.0  1.0  0.0 -1.0
%!       0.0  1.0 -2.0  0.0  3.0
%!       1.0  0.0  3.0 -1.0 -2.0];
%!
%! C = [ 1.0 -1.0  2.0 -2.0  0.0 -3.0
%!      -3.0  0.0  1.0 -1.0  1.0  0.0
%!       0.0  2.0  0.0 -4.0  0.0 -2.0
%!       1.0 -3.0  0.0  0.0  3.0  1.0
%!       0.0  1.0 -2.0  1.0  0.0 -2.0];
%!
%! D = [ 1.0 -1.0 -2.0  0.0  0.0
%!       0.0  1.0  0.0  1.0  0.0
%!       2.0 -1.0 -3.0  0.0  1.0
%!       0.0  1.0  0.0  1.0 -1.0
%!       0.0  0.0  1.0  2.0  1.0];
%!
%! P = ss (A, B, C, D, 1);
%! [K, ~, GAM] = hinfsyn (P, 2, 2, "gmax", 1000, "tolgam", 1e-4);
%! M = [K.A, K.B; K.C, K.D];
%!
%! KA = [-18.0030  52.0376  26.0831  -0.4271 -40.9022  18.0857
%!        18.8203 -57.6244 -29.0938   0.5870  45.3309 -19.8644
%!       -26.5994  77.9693  39.0368  -1.4020 -60.1129  26.6910
%!       -21.4163  62.1719  30.7507  -0.9201 -48.6221  21.8351
%!        -0.8911   4.2787   2.3286  -0.2424  -3.0376   1.2169
%!        -5.3286  16.1955   8.4824  -0.2489 -12.2348   5.1590];
%!
%! KB = [ 16.9788  14.1648
%!       -18.9215 -15.6726
%!        25.2046  21.2848
%!        20.1122  16.8322
%!         1.4104   1.2040
%!         5.3181   4.5149];
%!
%! KC = [ -9.1941  27.5165  13.7364  -0.3639 -21.5983   9.6025
%!         3.6490 -10.6194  -5.2772   0.2432   8.1108  -3.6293];
%!
%! KD = [  9.0317   7.5348
%!        -3.4006  -2.8219];
%!
%! M_exp = [KA, KB; KC, KD];
%! GAM_exp = 111.294;
%!
%!assert (M, M_exp, 1e-1);
%!assert (GAM, GAM_exp, 1e-3);
