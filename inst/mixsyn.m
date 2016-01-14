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
## @deftypefn{Function File} {[@var{K}, @var{N}, @var{gamma}, @var{info}] =} mixsyn (@var{G}, @var{W1}, @var{W2}, @var{W3}, @dots{})
## Solve stacked S/KS/T H-infinity problem.
## Mixed-sensitivity is the name given to transfer function shaping problems in which
## the sensitivity function
## @iftex
## @tex
## $ S = (I + G K)^{-1} $
## @end tex
## @end iftex
## @ifnottex
## @example
##              -1
## S = (I + G K)
## @end example
## @end ifnottex
## is shaped along with one or more other closed-loop transfer functions such as @var{K S}
## or the complementary sensitivity function
## @iftex
## @tex
## $ T = I - S = (I + G K)^{-1} G K $
## @end tex
## @end iftex
## @ifnottex
## @example
##                      -1
## T = I - S = (I + G K)
## @end example
## @end ifnottex
## in a typical one degree-of-freedom configuration, where @var{G} denotes the plant and
## @var{K} the (sub-)optimal controller to be found.  The shaping of multivariable
## transfer functions is based on the idea that a satisfactory definition of gain
## (range of gain) for a matrix transfer function is given by the singular values
## @iftex
## @tex
## $\\sigma$
## @end tex
## @end iftex
## @ifnottex
## 
## @end ifnottex
## of the transfer function. Hence the classical loop-shaping ideas of feedback design
## can be generalized to multivariable systems.  In addition to the requirement that
## @var{K} stabilizes @var{G}, the closed-loop objectives are as follows [1]:
## @enumerate
## @item For @emph{disturbance rejection} make
## @iftex
## @tex
## $\\overline{\\sigma}(S)$
## @end tex
## @end iftex
## @ifnottex
## 
## @end ifnottex
## small.
## @item For @emph{noise attenuation} make
## @iftex
## @tex
## $\\overline{\\sigma}(T)$
## @end tex
## @end iftex
## @ifnottex
## 
## @end ifnottex
## small.
## @item For @emph{reference tracking} make
## @iftex
## @tex
## $\\overline{\\sigma}(T) \\approx \\underline{\\sigma}(T) \\approx 1$.
## @end tex
## @end iftex
## @ifnottex
## 
## @end ifnottex
## @item For @emph{input usage (control energy) reduction} make
## @iftex
## @tex
## $\\overline{\\sigma}(K S)$
## @end tex
## @end iftex
## @ifnottex
## 
## @end ifnottex
## small.
## @item For @emph{robust stability} in the presence of an additive perturbation
## @iftex
## @tex
## $G_p = G + \\Delta$,
## @end tex
## @end iftex
## @ifnottex
## 
## @end ifnottex
## make
## @iftex
## @tex
## $\\overline{\\sigma}(K S)$
## @end tex
## @end iftex
## @ifnottex
## 
## @end ifnottex
## small.
## @item For @emph{robust stability} in the presence of a multiplicative output perturbation
## @iftex
## @tex
## $G_p = (I + \\Delta) G$,
## @end tex
## @end iftex
## @ifnottex
## 
## @end ifnottex
## make
## @iftex
## @tex
## $\\overline{\\sigma}(T)$
## @end tex
## @end iftex
## @ifnottex
## 
## @end ifnottex
## small.
## @end enumerate
## In order to find a robust controller for the so-called stacked
## @iftex
## @tex
## $S/KS/T \\ H_{\\infty}$
## @end tex
## @end iftex
## @ifnottex
## S/KS/T H-infinity
## @end ifnottex
## problem, the user function @command{mixsyn} minimizes the following criterion
## @iftex
## @tex
## $$ \\underset{K}{\\min} || N(K) ||_{\\infty}, \\quad N = | W_1 S; \\ W_2 K S; \\ W_3 T |$$
## @end tex
## @end iftex
## @ifnottex
## @example
##                              | W1 S   |
## min || N(K) ||           N = | W2 K S |
##  K            oo             | W3 T   |
## @end example
## @end ifnottex
## @code{[K, N] = mixsyn (G, W1, W2, W3)}.
## The user-defined weighting functions @var{W1}, @var{W2} and @var{W3} bound the largest
## singular values of the closed-loop transfer functions @var{S} (for performance),
## @var{K S} (to penalize large inputs) and @var{T} (for robustness and to avoid
## sensitivity to noise), respectively [1].
## A few points are to be considered when choosing the weights.
## The weigths @var{Wi} must all be proper and stable.  Therefore if one wishes,
## for example, to minimize @var{S} at low frequencies by a weighting @var{W1} including
## integral action, 
## @iftex
## @tex
## ${1 \\over s}$
## @end tex
## @end iftex
## @ifnottex
## @example
## 1
## -
## s
## @end example
## @end ifnottex
## needs to be approximated by
## @iftex
## @tex
## ${1 \\over s + \\epsilon}$, where $\\epsilon \\ll 1$.
## @end tex
## @end iftex
## @ifnottex
## @example
##   1
## -----    where   e << 1.
## s + e
## @end example
## @end ifnottex
## Similarly one might be interested in weighting @var{K S} with a non-proper weight
## @var{W2} to ensure that @var{K} is small outside the system bandwidth.
## The trick here is to replace a non-proper term such as
## @iftex
## @tex
## $ 1 + \\tau_1 s $ by $ {1 + \\tau_1 s \\over 1 + \\tau_2 s} $, where
## @end tex
## @end iftex
## @ifnottex
## @example
##                 1 + T1 s
## 1 + T1 s   by   --------,   where   T2 << T1.
##                 1 + T2 s
## @end example
## @end ifnottex
## @iftex
## @tex
## $\\tau_2 \\ll \\tau_1$
## @end tex
## @end iftex
## [1, 2]. 
##
##
## @strong{Inputs}
## @table @var
## @item G
## @acronym{LTI} model of plant.
## @item W1
## @acronym{LTI} model of performance weight.  Bounds the largest singular values of sensitivity @var{S}.
## Model must be empty @code{[]}, SISO or of appropriate size.
## @item W2
## @acronym{LTI} model to penalize large control inputs.  Bounds the largest singular values of @var{KS}.
## Model must be empty @code{[]}, SISO or of appropriate size.
## @item W3
## @acronym{LTI} model of robustness and noise sensitivity weight.  Bounds the largest singular values of 
## complementary sensitivity @var{T}.  Model must be empty @code{[]}, SISO or of appropriate size.
## @item @dots{}
## Optional arguments of @command{hinfsyn}.  Type @command{help hinfsyn} for more information.
## @end table
##
## All inputs must be proper/realizable.
## Scalars, vectors and matrices are possible instead of @acronym{LTI} models.
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
## @strong{Block Diagram}
## @example
## @group
##
##                                     | W1 S   |
## gamma = min||N(K)||             N = | W2 K S | = lft (P, K)
##          K         inf              | W3 T   |
## @end group
## @end example
## @example
## @group
##                                                       +------+  z1
##             +---------------------------------------->|  W1  |----->
##             |                                         +------+
##             |                                         +------+  z2
##             |                 +---------------------->|  W2  |----->
##             |                 |                       +------+
##  r   +    e |   +--------+  u |   +--------+  y       +------+  z3
## ----->(+)---+-->|  K(s)  |----+-->|  G(s)  |----+---->|  W3  |----->
##        ^ -      +--------+        +--------+    |     +------+
##        |                                        |
##        +----------------------------------------+
## @end group
## @end example
## @example
## @group
##                +--------+
##                |        |-----> z1 (p1x1)          z1 = W1 e
##  r (px1) ----->|  P(s)  |-----> z2 (p2x1)          z2 = W2 u
##                |        |-----> z3 (p3x1)          z3 = W3 y
##  u (mx1) ----->|        |-----> e (px1)            e = r - y
##                +--------+
## @end group
## @end example
## @example
## @group
##                +--------+  
##        r ----->|        |-----> z
##                |  P(s)  |
##        u +---->|        |-----+ e
##          |     +--------+     |
##          |                    |
##          |     +--------+     |
##          +-----|  K(s)  |<----+
##                +--------+
## @end group
## @end example
## @example
## @group
##                +--------+      
##        r ----->|  N(s)  |-----> z
##                +--------+
## @end group
## @end example
## @example
## @group
## Extended Plant:  P = augw (G, W1, W2, W3)
## Controller:      K = mixsyn (G, W1, W2, W3)
## Entire System:   N = lft (P, K)
## Open Loop:       L = G * K
## Closed Loop:     T = feedback (L)
## @end group
## @end example
##
## @strong{Algorithm}@*
## Relies on functions @command{augw} and @command{hinfsyn},
## which use SLICOT SB10FD, SB10DD and SB10AD by courtesy of
## @uref{http://www.slicot.org, NICONET e.V.}
##
## @strong{References}@*
## [1] Skogestad, S. and Postlethwaite I. (2005)
## @cite{Multivariable Feedback Control: Analysis and Design:
## Second Edition}.  Wiley, Chichester, England.@*
## [2] Meinsma, G. (1995)
## @cite{Unstable and nonproper weights in H-infinity control}
## Automatica, Vol. 31, No. 11, pp. 1655-1658
##
## @seealso{hinfsyn, augw}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: December 2009
## Version: 0.2

function [K, N, gamma, info] = mixsyn (G, W1 = [], W2 = [], W3 = [], varargin)

  if (nargin == 0)
    print_usage ();
  endif

  [p, m] = size (G);

  P = augw (G, W1, W2, W3);
  
  [K, N, gamma, info] = hinfsyn (P, p, m, varargin{:});

endfunction
