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
## @deftypefn {Function File} {[@var{scaledsys}, @var{info}] =} prescale (@var{sys})
## Prescale state-space model.
## Frequency response commands perform automatic scaling unless model property
## @var{scaled} is set to @var{true}.
##
## @strong{Inputs}
## @table @var
## @item sys
## LTI model.
## @end table
##
## @strong{Outputs}
## @table @var
## @item scaledsys
## Scaled state-space model.
## @item info
## Structure containing additional information.
## @item info.SL
## Left scaling factors.  @code{Tl = diag (info.SL)}.
## @item info.SR
## Right scaling factors.  @code{Tr = diag (info.SR)}.
## @end table
##
## @strong{Equations}
## @example
## @group
## Es = Tl * E * Tr
## As = Tl * A * Tr
## Bs = Tl * B
## Cs =      C * Tr
## Ds =      D
## @end group
## @end example
##
## For proper state-space models, @var{Tl} and @var{Tr} are inverse of each other.
##
## @strong{Algorithm}@*
## Uses SLICOT TB01ID and TG01AD by courtesy of
## @uref{http://www.slicot.org, NICONET e.V.}.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: June 2011
## Version: 0.1

function [retsys, varargout] = prescale (sys)

  if (nargin != 1)
    print_usage ();
  endif

  if (! isa (sys, "ss"))
    warning ("prescale: system not in state-space form");
    sys = ss (sys);
  endif

  [retsys, lscale, rscale] = __prescale__ (sys);
  
  if (nargout > 1)
    varargout{1} = struct ("SL", lscale, "SR", rscale);
  endif

endfunction
