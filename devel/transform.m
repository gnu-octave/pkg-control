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
## @deftypefn {Function File} {@var{sys} =} transform (@var{sys}, @var{Tx})
## @deftypefnx {Function File} {@var{sys} =} transform (@var{sys}, @var{Tx}, @var{Tu}, @var{Ty})
## @deftypefnx {Function File} {@var{sys} =} transform (@var{sys}, @var{Tl}, @var{Tr})
## @deftypefnx {Function File} {@var{sys} =} transform (@var{sys}, @var{Tl}, @var{Tr}, @var{Tu}, @var{Ty})
## Transform state-space model.
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
## x = Tx * xt
## u = Tu * ut
## y = Ty * yt
##
## At = Tx \ A * Tx
## Bt = Tx \ B * Tu
## Ct = Ty \ C * Tx
## Dt = Ty \ D * Tu
##
## Et = Tl * E * Tr
## At = Tl * A * Tr
## Bt = Tl * B * Tu
## Ct = Ty \ C * Tr
## Dt = Ty \ D * Tu
##
## For proper state-space models, Tl and Tr are inverse of each other.
## @end group
## @end example
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2011
## Version: 0.1

## To be placed inside @lti

function retsys = transform (sys, Tl, Tr = [], Tu = [], Ty = [])

  if (nargin < 2 || nargin > 5)
    print_usage ();
  endif
  
  if (! isa (sys, "ss"))
    warning ("transform: system not in state-space form");
    sys = ss (sys);
  endif
  
  if (! is_real_square_matrix (Tl, Tr, Tu, Ty))
    error ("transform: transformation matrices must be square and real");
  endif

  [A, B, C, D, E, tsam] = dssdata (sys, []);

  if (isempty (E))
    Tx = Tl;
    Ty = Tu;
    Tu = Tr;
  endif

  if (isempty (Tu))
    Tu = eye (columns (B));
  endif

  if (isempty (Ty))
    Ty = eye (rows (C));
  endif

  if (isempty (E))
    At = Tx \ A * Tx;
    Bt = Tx \ B * Tu;
    Ct = Ty \ C * Tx;
    Dt = Ty \ D * Tu;
    retsys = ss (At, Bt, Ct, Dt, tsam);
  else
    Et = Tl * E * Tr;
    At = Tl * A * Tr;
    Bt = Tl * B * Tu;
    Ct = Ty \ C * Tr;
    Dt = Ty \ D * Tu;
    retsys = dss (At, Bt, Ct, Dt, Et, tsam);
  endif

endfunction









