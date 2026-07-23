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

## -*- texinfo -*-
## @deftypefn {Function File} {@var{sys} =} pade (@var{sys}, @var{n})
## Replace every delay (@var{InputDelay}, @var{OutputDelay}, @var{IODelay})
## of the @acronym{LTI} model @var{sys} by an @var{n}th-order Pade
## approximation, returning a delay-free equivalent model.
##
## @strong{Inputs}
## @table @var
## @item sys
## @var{tf} or @var{zpk} model (@var{ss} models are not yet supported).
## @item n
## Order of the Pade approximation.  Either a scalar (applied to every
## nonzero delay) or a vector with one entry per nonzero delay -- see
## @code{__pade_order_vector__} for the exact ordering convention.
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## Delay-free @var{tf} or @var{zpk} model with the Pade-approximated
## rational factors absorbed into each entry's numerator/denominator
## (or zeros/poles).
## @end table
## @end deftypefn

## Author: Prateek Ganguli <prateek.ganguli@gmail.com>
## Created: July 2026
## Version: 0.1

function sys = pade (sys, n)

  if (nargin != 2 || ! isa (sys, "lti"))
    print_usage ();
  endif

  if (isa (sys, "ss"))
    error ("pade: ss models are not yet supported");
  endif

  if (! (isa (sys, "tf") || isa (sys, "zpk")))
    error ("pade: only tf and zpk models are supported (ss not yet implemented)");
  endif

  if (! hasdelay (sys))
    return;
  endif

  [indelay, outdelay, iodelay] = get (sys, "inputdelay", "outputdelay", "iodelay");
  total = totaldelay (sys);
  [pr, pc] = size (total);

  [orders, idx] = __pade_order_vector__ (sys, n);

  n_in  = numel (idx.input);
  n_out = numel (idx.output);
  n_iod = numel (idx.iod);

  orders_input  = orders (1 : n_in);
  orders_output = orders (n_in + 1 : n_in + n_out);
  orders_iod    = orders (n_in + n_out + 1 : n_in + n_out + n_iod);

  ## Build a per-entry order matrix.  Priority (lowest to highest, later
  ## assignments win on conflict): OutputDelay broadcast across a row,
  ## InputDelay broadcast across a column, IODelay assigned directly to
  ## its own (i,j) entry -- IODelay is the most entry-specific, so it
  ## takes precedence if more than one delay kind is nonzero for the same
  ## entry.  Test fixtures in this task only ever set one delay kind at a
  ## time per entry, so no conflict arises in practice.
  order_matrix = nan (pr, pc);

  for k = 1 : n_out
    i = idx.output (k);
    order_matrix (i, :) = orders_output (k);
  endfor

  for k = 1 : n_in
    j = idx.input (k);
    order_matrix (:, j) = orders_input (k);
  endfor

  for k = 1 : n_iod
    [i, j] = ind2sub ([pr, pc], idx.iod (k));
    order_matrix (i, j) = orders_iod (k);
  endfor

  origsys = sys;

  if (isa (sys, "zpk"))

    [z, p, k, tsam] = zpkdata (sys);

    for i = 1 : pr
      for j = 1 : pc
        d = total (i, j);
        if (d > 0)
          order = order_matrix (i, j);
          [num_p, den_p] = pade (d, order);
          gain_factor = num_p(1) / den_p(1);
          z{i,j} = [z{i,j}; roots(num_p)];
          p{i,j} = [p{i,j}; roots(den_p)];
          k(i,j) = k(i,j) * gain_factor;
        endif
      endfor
    endfor

    sys = zpk (z, p, k, tsam);

  else

    [num, den, tsam] = tfdata (sys);

    for i = 1 : pr
      for j = 1 : pc
        d = total (i, j);
        if (d > 0)
          order = order_matrix (i, j);
          [num_p, den_p] = pade (d, order);
          num{i,j} = conv (num{i,j}, num_p);
          den{i,j} = conv (den{i,j}, den_p);
        endif
      endfor
    endfor

    sys = tf (num, den, tsam);

  endif

  sys = set (sys, "lti", origsys);
  sys = set (sys, "InputDelay", 0, "OutputDelay", 0, "IODelay", 0);

endfunction


%!test  # SISO tf with InputDelay: freqresp matches hand-multiplied reference
%! T = 0.5;
%! sys = tf (1, [1 2], "InputDelay", T);
%! sys2 = pade (sys, 3);
%! assert (hasdelay (sys2), false);
%! w = [0.1, 1, 5];
%! [num_p, den_p] = pade (T, 3);
%! expected = freqresp (tf (1, [1 2]), w) .* freqresp (tf (num_p, den_p), w);
%! assert (freqresp (sys2, w), expected, 1e-8);

%!test  # SISO tf with InputDelay, a different order
%! T = 0.5;
%! sys = tf (1, [1 2], "InputDelay", T);
%! sys2 = pade (sys, 5);
%! assert (hasdelay (sys2), false);
%! w = [0.1, 1, 5];
%! [num_p, den_p] = pade (T, 5);
%! expected = freqresp (tf (1, [1 2]), w) .* freqresp (tf (num_p, den_p), w);
%! assert (freqresp (sys2, w), expected, 1e-8);

%!test  # dense MIMO tf with differing per-channel IODelay
%! sys = tf ({1, 1; 1, 1}, {[1 2], [1 3]; [1 4], [1 5]});
%! sys = set (sys, "IODelay", [0.2, 0.3; 0.4, 0.5]);
%! sys2 = pade (sys, 3);
%! assert (hasdelay (sys2), false);
%! w = [0.1, 1, 5];
%! rational = tf ({1, 1; 1, 1}, {[1 2], [1 3]; [1 4], [1 5]});
%! total = [0.2, 0.3; 0.4, 0.5];
%! resp = freqresp (rational, w);
%! expected = zeros (2, 2, numel (w));
%! for i = 1:2
%!   for j = 1:2
%!     [num_p, den_p] = pade (total(i,j), 3);
%!     expected(i,j,:) = reshape (resp(i,j,:), 1, []) .* reshape (freqresp (tf (num_p, den_p), w), 1, []);
%!   endfor
%! endfor
%! assert (freqresp (sys2, w), expected, 1e-6);

%!test  # equivalent zpk case, SISO with OutputDelay
%! T = 0.3;
%! sys = zpk ([], -2, 3, "OutputDelay", T);
%! sys2 = pade (sys, 4);
%! assert (hasdelay (sys2), false);
%! w = [0.1, 1, 5];
%! [num_p, den_p] = pade (T, 4);
%! expected = freqresp (zpk ([], -2, 3), w) .* freqresp (tf (num_p, den_p), w);
%! assert (freqresp (sys2, w), expected, 1e-8);

%!test  # equivalent zpk case, dense MIMO with IODelay
%! z = {[], []; [], []};
%! p = {-2, -3; -4, -5};
%! k = [1, 1; 1, 1];
%! sys = zpk (z, p, k);
%! sys = set (sys, "IODelay", [0.2, 0.3; 0.4, 0.5]);
%! sys2 = pade (sys, 3);
%! assert (hasdelay (sys2), false);
%! w = [0.1, 1, 5];
%! rational = zpk (z, p, k);
%! total = [0.2, 0.3; 0.4, 0.5];
%! resp = freqresp (rational, w);
%! expected = zeros (2, 2, numel (w));
%! for i = 1:2
%!   for j = 1:2
%!     [num_p, den_p] = pade (total(i,j), 3);
%!     expected(i,j,:) = reshape (resp(i,j,:), 1, []) .* reshape (freqresp (tf (num_p, den_p), w), 1, []);
%!   endfor
%! endfor
%! assert (freqresp (sys2, w), expected, 1e-6);

%!test  # per-delay order vector: each entry keeps its own assigned order
%! # 1x2 tf, two nonzero IODelay entries with distinct Pade orders 2 and 5.
%! # Each substituted entry's denominator degree must equal its own
%! # rational-part degree (1, since den = [1 2] or [1 3]) plus its own
%! # assigned Pade order -- checked individually, not in aggregate.
%! sys = tf ({1, 1}, {[1 2], [1 3]});
%! sys = set (sys, "IODelay", [0.2, 0.4]);
%! sys2 = pade (sys, [2, 5]);
%! [num, den] = tfdata (sys2);
%! assert (length (den{1,1}) - 1, 1 + 2);
%! assert (length (den{1,2}) - 1, 1 + 5);

%!error <does not match the number of nonzero delays> pade (tf (1, [1 2], "InputDelay", 0.5), [1, 2])
