## Copyright (C) 2012   Lukas F. Reichlin
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
## @deftypefn{Function File} {@var{Ms} =} sensitivity (@var{L})
## @deftypefnx{Function File} {@var{Ms} =} sensitivity (@var{P}, @var{C})
## @deftypefnx{Function File} {@var{Ms} =} sensitivity (@var{P}, @var{C1}, @var{C2}, @dots{})
## Return sensitivity margin @var{Ms}.
##
## @strong{Inputs}
## @table @var
## @item L
## Open loop transfer function.
## @var{L} can be any type of LTI system, but it must be square.
## @item P
## Plant model.  Any type of LTI system.
## @item C
## Controller model.  Any type of LTI system.
## @item C1, C2, @dots{}
## If several controllers are specified, command @command{sensitivity}
## computes the sensitivity @var{Ms} for each of them in combination
## with plant @var{P}.
## @end table
##
## @strong{Outputs}
## @table @var
## @item Ms
## Sensitivity margin @var{Ms} as defined in [1].  Scalar value.
## If several controllers are specified, @var{Ms} becomes
## a vector with as many entries as controllers.
## @end table
##
## @strong{Algorithm}@*
## Uses SLICOT AB13DD by courtesy of
## @uref{http://www.slicot.org, NICONET e.V.}
##
## @strong{References}@*
## [1] Astr@"om, K. and H@"agglund, T. (1995)
## PID Controllers:
## Theory, Design and Tuning,
## Second Edition.
## Instrument Society of America.
##
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: August 2012
## Version: 0.1

function Ms = sensitivity (G, varargin)

  ## TODO: show nyquist diagram of L with circle
  ##       center (-1, 0) and radius equal to the
  ##       shortest distance 1/Ms, frequency w is
  ##       [Ms, w] = norm (S, inf)

  if (nargin == 0)
    print_usage ();
  elseif (nargin == 1)              # L := G
    I = eye (size (G));
    S = feedback (I, G);            # S = inv (I + G),  S = feedback (I, G*-I, "+")
    Ms = norm (S, inf);
  else                              # P := G,  C := varargin
    L = cellfun (@(C) G*C, varargin, "uniformoutput", false);
    I = cellfun (@(L) eye (size (L)), L, "uniformoutput", false);
    S = cellfun (@feedback, I, L, "uniformoutput", false);
    Ms = cellfun (@(S) norm (S, inf), S);
  endif

endfunction
