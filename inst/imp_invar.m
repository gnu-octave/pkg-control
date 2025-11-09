## Copyright (C) 2009-2015   Lukas F. Reichlin
## Copyright (C) 2016 Douglas A. Stewart
## Copyright (C) 2024 Torsten Lilge
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
## along with LTI Syncope. If not, see <http://www.gnu.org/licenses/>.


## -*- texinfo -*-
## @deftypefn  {Function File} {[@var{b_out}, @var{a_out}] =} imp_invar (@var{b}, @var{a}, @var{fs}, @var{tol})
## @deftypefnx {Function File} {[@var{b_out}, @var{a_out}] =} imp_invar (@var{b}, @var{a}, @var{fs})
## @deftypefnx {Function File} {[@var{b_out}, @var{a_out}] =} imp_invar (@var{b}, @var{a})
## @deftypefnx {Function File} {[@var{sys_out}] =} imp_invar (@var{b}, @var{a}, @var{fs}. @var{tol})
## @deftypefnx {Function File} {[@var{sys_out}] =} imp_invar (@var{sys_in}, @var{fs}, @var{tol})
## Converts analog filter with coefficients @var{b} and @var{a} and/or @var{sys_in} to digital,
## conserving impulse response.
##
## MIMO systems are only supported with @var{sys_in} as input argument.
##
## @strong{Inputs}
## @table @var
## @item b
## Numerator coefficients of continuous-time LTI system. 
## @item a
## Denominator coefficients of continuous-time LTI system.
## @item fs 
## Sampling frequency. If @var{fs} is not specified, or is an empty vector,
## it defaults to 1Hz.
## tol
## Tolerance of the internally used function minreal for canceling identical
## poles and zeros. If @var{tol} is not specified, it defaults to 0.0001 (0.1%).
## @item sys_in
## System definition of the continuous-time LTI system. This can also be
## a MIMO system.
## @end table
##
## @strong{Outputs}
## @table @var
## @item b_out
## Numerator coefficients of the discrete-time impulse invariant LTI system. 
## @item a_out
## Denominator coefficients of the discrete-time impulse invariant LTI system. 
## @item sys_out
## Discrete-time impulse invaraiant LTI system. If @var{sys_in} is given as
## state space representation, @var{sys_out} is also returned in state space,
## otherwise as transfer function.
## @end table
##
## @strong{Algorithm}
##
## The step equivalent discretization of G(s) (zoh) results in
## G_zoh(z) = (z-1)/z * Z@{G(s)/s@} where Z@{@} is the z-transformation.
## The transfer function of the impulse equivalent discretization
## is given by T*Z@{G(s)@}. Therefore, the zoh discretizaiton method for
## s*G(s) multipled by T*z/(z-1) leads to the desired result.
##
## @strong{Remark}
##
## For the impulse response of a discrete-time system, the input
## sequence @{1/T,0,0,0,...@} and not the unit impulse is considered.
## For this reason, the factor T is required for the impulse invaraint
## discretizaiton (see Algorithm).
##
## @seealso{@@lti/c2d}
## @end deftypefn


function [bz az] = imp_invar (b , a , fs , tol = 1e-4)
  ## This funtion will accept both a
  ## sys variable as input and/or
  ## numerator, denominator as input.

  if (nargin < 1)
    print_usage;
  endif

  if (isa (b, "lti"))

    ## the input is an LTI object, therefore inputs are (sys,fs,tol)
    ## so b is sys, a is fs, and fs is tol
    ## in this case, MIMO systems are allowed
    if (exist("fs","var") != 0)
      tol = fs;
    else
      tol=0.0001;
    endif

    if (exist ("a") == 1)
      fs=a;
    else
      fs=1;
    endif

    ## lti system given in state space and nargout == 1?
    ## if yes, just return discretized system in state space
    if ( isa (b,"ss") && (nargout == 1) )
      bz = c2d (b, 1/fs, "imp");
      return;
    endif

    ## get all polynomials in a cell array
    [ny, nu] = size (b);
    [bcell, acell] = tfdata (b);

  else

    ## some internal functions call imp_invar with polynomials in cells
    if (iscell (b))
      b = b{1,1};
    endif
    if (iscell (a))
      a = a{1,1};
    endif

    ## the input is vectors
    if (! (ismatrix (b) && ismatrix (a))) || ...
        ((min (size (b)) != 1) && (min (size (a)) != 1))
      error ("imp_invar: first two arguments must be vectors\n");
    endif
    if (exist ("fs") == 0)
      fs = 1;
    endif

    ny = nu = 1;
    bcell = cell ();
    acell = cell ();
    bcell{1,1} = b;
    acell{1,1} = a;

  endif

  if (isempty (fs))
    fs = 1;
  endif

  if (isempty (tol))
    tol = 1e-4;
  endif

  T = 1/fs;

  bz = cell (ny,nu);
  az = bz;

  for iy = 1:ny
    for iu = 1:nu

      b = remove_leading_zeros (bcell{iy,iu});
      a = remove_leading_zeros (acell{iy,iu});

      if (length (b) >= length (a))
        error("Order numerator >= order denominator");
      endif

      ## Apply zoh method for s*G(s) and multiply the result by z/(z-1).
      b = conv (b, [1 0]);              # multiply by s
      G_zoh = c2d (tf (b,a), T, 'zoh'); # zoh method for s*G(s)

      [bzz,azz] = tfdata (G_zoh, 'v');  # get polynomials of result
      bzz = remove_leading_zeros (bzz);
      bzz = conv (bzz, [T 0]);          # multiply numerator by T*z
      azz = conv (azz, [1 -1]);         # multiply denominator by z-1

      sys1 = tf (bzz, azz, T);
      sys2 = minreal (sys1, tol); # Use this to remove the common roots.

      [bz{iy,iu}, az{iy,iu}] = tfdata (sys2, "v");

    endfor
  endfor

  if (nargout() < 2)
    bz = tf (bz, az, T);
  else
    if (ny*nu == 1)
      bz = bz{1,1};
      az = az{1,1};
    endif
  endif

endfunction


function x_clean = remove_leading_zeros (x)

  nonzero = find (x);
  if length (nonzero) == 0
    x_clean = 0;
  else
    x_clean = x(nonzero(1):end);
  endif

endfunction


## Tests
##
%!shared bz1, az1, bz2, az2, bz1_e, az1_e, bz2_e, az2_e
%!
%! s = tf ('s');
%! Gs = (s-2)*(s-1)*(s+5)/s/(s+1)/(s+2)^3/(s+3)/(s+4);
%! [b,a] = tfdata (Gs, 'v');
%! [bz1,az1] = imp_invar (Gs, 2);
%! [bz2,az2] = imp_invar (b, a, 5);
%!
%! bz1_e = 1/2*[-0.0000  0.0036 -0.0128  0.0039  0.0125 -0.0001 -0.0001  0.0000];
%! az1_e = [ 1.0000 -3.0686  3.7873 -2.4518  0.9020 -0.1886  0.0207 -0.0009];
%!
%! bz2_e = 1/5*[-0.0000  0.0007 -0.0007  -0.0025  0.0032 -0.0004 -0.0001  0.0000];
%! az2_e = [ 1.0000 -4.8278  9.8933 -11.1569  7.4787 -2.9798  0.6534 -0.0608];
%!
%!assert (az1, az1_e, 1e-4);
%!assert (bz1, bz1_e, 1e-4);
%!assert (az2, az2_e, 1e-4);
%!assert (bz2, bz2_e, 1e-4);

