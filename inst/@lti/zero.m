## Copyright (C) 2009, 2011   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{z} =} zero (@var{sys})
## @deftypefnx {Function File} {[@var{z}, @var{k}] =} zero (@var{sys})
## Compute invariant zeros and gain of LTI model.
## Invariant zeros are also known as Smith zeros.
##
## @strong{Inputs}
## @table @var
## @item sys
## LTI model.
## @end table
##
## @strong{Outputs}
## @table @var
## @item z
## Invariant zeros of @var{sys} as defined in [1].
## @item k
## Gain of @var{sys}.
## @end table
##
## @strong{Algorithm}@*
## For (descriptor) state-space models, @command{zero}
## relies on SLICOT AB08ND and AG08BD by courtesy of
## @uref{http://www.slicot.org, NICONET e.V.}
## For @acro{SISO} transfer functions, @command{zero}
## uses Octave's @command{roots}.
## @acro{MIMO} transfer functions are converted to
## a minimal state-space representation for the
## computation of the zeros.
##
## @strong{Reference}@*
## [1] MacFarlane, A. and Karcanias, N.
## @cite{Poles and zeros of linear multivariable systems:
## a survey of the algebraic, geometric and complex-variable
## theory}.  Int. J. Control, vol. 24, pp. 33-74, 1976.@*
##
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.3

function [zer, gain] = zero (sys)

  if (nargin > 1)
    print_usage ();
  endif

  [zer, gain] = __zero__ (sys, nargout);

endfunction


## Invariant zeros of state-space models
##
## Results from the "Dark Side" 7.5 and 7.8
##
##  -13.2759
##   12.5774
##  -0.0155
##
## Results from Scilab 5.2.0b1 (trzeros)
##
##  - 13.275931  
##    12.577369  
##  - 0.0155265
##
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
%! sys = ss (A, B, C, D, "scaled", true);
%! z = sort (zero (sys));
%!
%! z_exp = sort ([-13.2759; 12.5774; -0.0155]);
%!
%!assert (z, z_exp, 1e-4);


## Invariant zeros of descriptor state-space models
%!shared z, z_exp
%! A = [  1     0     0     0     0     0     0     0     0
%!        0     1     0     0     0     0     0     0     0
%!        0     0     1     0     0     0     0     0     0
%!        0     0     0     1     0     0     0     0     0
%!        0     0     0     0     1     0     0     0     0
%!        0     0     0     0     0     1     0     0     0
%!        0     0     0     0     0     0     1     0     0
%!        0     0     0     0     0     0     0     1     0
%!        0     0     0     0     0     0     0     0     1 ];
%!
%! E = [  0     0     0     0     0     0     0     0     0
%!        1     0     0     0     0     0     0     0     0
%!        0     1     0     0     0     0     0     0     0
%!        0     0     0     0     0     0     0     0     0
%!        0     0     0     1     0     0     0     0     0
%!        0     0     0     0     1     0     0     0     0
%!        0     0     0     0     0     0     0     0     0
%!        0     0     0     0     0     0     1     0     0
%!        0     0     0     0     0     0     0     1     0 ];
%!
%! B = [ -1     0     0
%!        0     0     0
%!        0     0     0
%!        0    -1     0
%!        0     0     0
%!        0     0     0
%!        0     0    -1
%!        0     0     0
%!        0     0     0 ];
%!
%! C = [  0     1     1     0     3     4     0     0     2
%!        0     1     0     0     4     0     0     2     0
%!        0     0     1     0    -1     4     0    -2     2 ];
%!
%! D = [  1     2    -2
%!        0    -1    -2
%!        0     0     0 ];
%!
%! sys = dss (A, B, C, D, E, "scaled", true);
%! z = zero (sys);
%!
%! z_exp = 1;
%!
%!assert (z, z_exp, 1e-4);


## Gain of descriptor state-space models
%!shared p, pi, z, zi, k, ki, p_tf, pi_tf, z_tf, zi_tf, k_tf, ki_tf
%! P = ss (-2, 3, 4, 5);
%! Pi = inv (P);
%!
%! p = pole (P);
%! [z, k] = zero (P);
%!
%! pi = pole (Pi);
%! [zi, ki] = zero (Pi);
%!
%! P_tf = tf (P);
%! Pi_tf = tf (Pi);
%!
%! p_tf = pole (P_tf);
%! [z_tf, k_tf] = zero (P_tf);
%!
%! pi_tf = pole (Pi_tf);
%! [zi_tf, ki_tf] = zero (Pi_tf);
%!
%!assert (p, zi, 1e-4);
%!assert (z, pi, 1e-4);
%!assert (k, inv (ki), 1e-4);
%!assert (p_tf, zi_tf, 1e-4);
%!assert (z_tf, pi_tf, 1e-4);
%!assert (k_tf, inv (ki_tf), 1e-4);


## Example taken from Paper [1]
%!shared zo, ze
%! A = diag ([1, 1, 3, -4, -1, 3]);
%! 
%! B = [  0,  -1
%!       -1,   0
%!        1,  -1
%!        0,   0
%!        0,   1
%!       -1,  -1  ];
%!        
%! C = [  1,  0,  0,  1,  0,  0
%!        0,  1,  0,  1,  0,  1
%!        0,  0,  1,  0,  0,  1  ];
%!         
%! D = zeros (3, 2);
%! 
%! SYS = ss (A, B, C, D);
%! zo = zero (SYS);
%! 
%! ze = [2; -1];
%! 
%!assert (zo, ze, 1e-4); 
