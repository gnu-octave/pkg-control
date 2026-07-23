## Copyright (C) 2026        Prateek Ganguli
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

## Compute the per-nonzero-delay Pade approximation order vector for
## pade(sys, n).
##
## Shared helper for @lti/pade.m (tf/zpk) and the future ss-specific pade
## dispatches (Tasks 2/3) -- build the per-entry order list once here and
## have every case index into it consistently.
##
## Ordering convention (this codebase's own, since MATLAB's internal
## ordering isn't documented/verifiable): the total nonzero-delay count,
## and the order in which "n" (if a vector) is consumed, walks:
##   1. InputDelay nonzero entries, in input order.
##   2. OutputDelay nonzero entries, in output order.
##   3. IODelay nonzero entries, in column-major (Octave linear index) order.
##   4. InternalDelay ports, in the existing internal port order returned
##      by get(sys, "internaldelay") (always empty for tf/zpk, since those
##      types can never carry InternalDelay).
##
## Inputs:
##   sys - an LTI model.
##   n   - scalar or vector Pade order.  If scalar, the same order is
##         used for every nonzero delay.  If a vector, its length must
##         equal the total nonzero-delay count (computed per the ordering
##         above); anything else is an error.
##
## Outputs:
##   orders - column vector of per-delay-entry orders, length equal to the
##            total nonzero-delay count, in the documented order above.
##   idx    - struct with fields "input", "output", "iod", "internal", each
##            a vector of linear/element indices (into InputDelay,
##            OutputDelay, IODelay(:), and internaldelay(:) respectively)
##            identifying which entries are nonzero, in the same order
##            "orders" was built in -- so callers can walk
##            idx.input/idx.output/idx.iod/idx.internal in lockstep with
##            orders(1:numel(idx.input)), orders(numel(idx.input)+1 : ...),
##            etc.
function [orders, idx] = __pade_order_vector__ (sys, n)

  if (nargin != 2 || ! isa (sys, "lti"))
    print_usage ();
  endif

  [indelay, outdelay, iodelay] = get (sys, "inputdelay", "outputdelay", "iodelay");

  ## InternalDelay is a ss-specific property; get (sys, "internaldelay")
  ## errors for tf/zpk models, so only query it for ss (always empty for
  ## tf/zpk, matching the documented ordering convention above).
  if (isa (sys, "ss"))
    internaldelay = get (sys, "internaldelay");
  else
    internaldelay = [];
  endif

  idx.input    = find (indelay(:)    != 0);
  idx.output   = find (outdelay(:)   != 0);
  idx.iod      = find (iodelay(:)    != 0);
  idx.internal = find (internaldelay(:) != 0);

  total = numel (idx.input) + numel (idx.output) + numel (idx.iod) + numel (idx.internal);

  if (isscalar (n))
    orders = repmat (n, total, 1);
  else
    if (numel (n) != total)
      error ("pade: order vector length (%d) does not match the number of nonzero delays (%d)",
             numel (n), total);
    endif
    orders = n(:);
  endif

endfunction
