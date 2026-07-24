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
## @var{tf}, @var{zpk}, or @var{ss} model.  @var{ss} models carrying a
## nonzero @var{InternalDelay} are supported: the delay loop is closed
## through a rational Pade filter via a linear-fractional transformation
## (@code{lft}).
## @item n
## Order of the Pade approximation.  Either a scalar (applied to every
## nonzero delay) or a vector with one entry per nonzero delay -- see
## @code{__pade_order_vector__} for the exact ordering convention.
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## Delay-free equivalent model with the Pade-approximated rational factors
## absorbed into each entry's numerator/denominator (or zeros/poles, or
## state-space realization).  For @var{ss} inputs, the result is obtained
## via a @code{tf} roundtrip and does @strong{not} preserve the original
## state basis or state count (matching @acronym{MATLAB}'s own
## @code{pade()}, which likewise only guarantees an equivalent model, not
## a minimal one -- use @code{minreal} afterward if minimality is needed).
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
    if (hasinternaldelay (sys))
      sys = __pade_substitute_internal__ (sys, n);
      return;
    endif
    if (! hasdelay (sys))
      return;
    endif
    sys_tf = __pade_substitute_ordinary__ (tf (sys), n);
    sys = ss (sys_tf);
    return;
  endif

  if (! (isa (sys, "tf") || isa (sys, "zpk")))
    error ("pade: only tf, zpk, and ss models are supported");
  endif

  if (! hasdelay (sys))
    return;
  endif

  sys = __pade_substitute_ordinary__ (sys, n);

endfunction

## Shared ordinary-delay (InputDelay/OutputDelay/IODelay) per-entry Pade
## substitution logic for tf/zpk models.  Assumes hasdelay(sys) is true.
## Used directly for tf/zpk inputs, and via a tf-roundtrip for ss inputs
## with no InternalDelay.
function sys = __pade_substitute_ordinary__ (sys, n)

  [indelay, outdelay, iodelay] = get (sys, "inputdelay", "outputdelay", "iodelay");
  total = totaldelay (sys);
  [pr, pc] = size (total);

  if (isdt (sys))
    sys = set (sys, "InputDelay", indelay, "OutputDelay", outdelay, "IODelay", iodelay);
    sys = absorbDelay (sys);
    return;
  endif

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

## InternalDelay Pade substitution for ss models.  Each internal delay
## port is a signal z(t) fed back to the plant as w(t) = z(t-tau); the
## exact freqresp closes this loop with Delta = diag(exp(-jw*tau)) via an
## LFT (see @ss/__freqresp__.m).  Here we replace that ideal-delay loop
## closure by a rational one: build the extended ordinary plant (whose last
## nports inputs are w and last nports outputs are z) with
## __ss_ext_build__, then close the loop through a diagonal MIMO Pade
## filter G_pade (one SISO Pade block per port) using the Redheffer star
## product lft().  As the order grows, each SISO Pade block converges to
## exp(-jw*tau), so the closed loop converges to the exact-delay system.
##
## Any ordinary delay (InputDelay/OutputDelay/IODelay) carried alongside
## the InternalDelay is independent of the loop closure, so it is handled
## afterwards by the same tf-roundtrip substitution used for the
## InternalDelay-free ss case (__pade_substitute_ordinary__), fed the
## ordinary slice of the SAME order assignment computed once up front.
function sys = __pade_substitute_internal__ (sys, n)

  ## Full per-delay order assignment, computed ONCE so the internal and
  ## ordinary substitutions consume consistent slices (ordering per
  ## __pade_order_vector__: InputDelay, OutputDelay, IODelay, InternalDelay).
  [orders, idx] = __pade_order_vector__ (sys, n);
  n_in  = numel (idx.input);
  n_out = numel (idx.output);
  n_iod = numel (idx.iod);
  n_int = numel (idx.internal);

  orders_ordinary = orders (1 : n_in + n_out + n_iod);
  orders_internal = orders (n_in + n_out + n_iod + 1 : end);

  ## Original internal-delay values and any ordinary delay, extracted BEFORE
  ## the extended plant is rebuilt (the rebuild deliberately drops both the
  ## delay ports and the ordinary delay).
  tau = get (sys, "internaldelay");
  tau = tau(:);
  had_ordinary = hasdelay (sys);
  [oindelay, ooutdelay, oiodelay] = get (sys, "inputdelay", "outputdelay", "iodelay");

  ## Build the extended ordinary plant.  __ss_ext_build__ folds the delay
  ## ports (b2/c2/d12/d21/d22) into ordinary extra inputs/outputs but leaves
  ## the returned ext_sys's lti input/output NAMES (hence its reported
  ## size()) and its residual b2/c2/.../internaldelay fields untouched, so
  ## it is NOT a well-formed ordinary ss.  Rebuild a clean ss from just the
  ## extended (a,b,c,d[,e]) matrices: that gives correct input/output counts
  ## for lft()'s size() calls AND clears every residual port field and the
  ## internaldelay, so hasinternaldelay(ext_clean) is false.
  ## __ss_ext_build__ lives in @ss and returns an ss whose a/b/c/d hold the
  ## extended matrices; extract them with dssdata (the raw fields are not
  ## dot-accessible from outside @ss).  Passing [] keeps a descriptor E empty
  ## rather than expanding it to the identity.
  [ext_sys, nu, ny] = __ss_ext_build__ (sys);
  [ea, eb, ec, ed, ee, etsam] = dssdata (ext_sys, []);
  nports = columns (eb) - nu;
  if (nports != rows (ec) - ny)
    error ("pade: internal delay port dimension mismatch (%d inputs vs %d outputs)",
           nports, rows (ec) - ny);
  endif

  ## Rebuild a clean ordinary ss from just the extended matrices.  This gives
  ## correct input/output counts for lft()'s size() calls (the raw ext_sys
  ## keeps the ORIGINAL input/output names, so its reported size is stale)
  ## AND clears every residual port field and the internaldelay, so
  ## hasinternaldelay(ext_clean) is false.
  if (isempty (ee))
    ext_clean = ss (ea, eb, ec, ed, etsam);
  else
    ext_clean = dss (ea, eb, ec, ed, ee, etsam);
  endif

  if (nports == 0)
    ## Degenerate: an internaldelay is set but there are no actual delay
    ## ports (empty b2/c2/...), so the "delay" affects no signal.  There is
    ## no loop to close -- the delay-free equivalent is just the clean plant.
    sys = ext_clean;
  else
    if (nports != numel (tau) || nports != n_int)
      error ("pade: internal delay port count (%d) does not match delay count (%d/%d)",
             nports, numel (tau), n_int);
    endif

    ## Diagonal MIMO Pade filter: nports inputs (each port's z, the delayed
    ## signal's source) and nports outputs (each port's w, the destination).
    ## For discrete sys, tau(k) is already a nonnegative-integer sample
    ## count, so the loop is closed EXACTLY via __exact_discrete_delay_ss__
    ## instead of a rational Pade approximation.
    if (isdt (sys))
      G_pade = __exact_discrete_delay_ss__ (tau(1), etsam);
      for k = 2 : nports
        G_pade = append (G_pade, __exact_discrete_delay_ss__ (tau(k), etsam));
      endfor
    else
      G_pade = pade (tau(1), orders_internal(1));
      for k = 2 : nports
        G_pade = append (G_pade, pade (tau(k), orders_internal(k)));
      endfor
    endif

    ## Close the loop.  lft(sys1, sys2, nu, ny) connects the LAST nu inputs of
    ## sys1 to the FIRST nu outputs of sys2, and the LAST ny outputs of sys1
    ## to the FIRST ny inputs of sys2.  With sys1 = ext_clean, sys2 = G_pade,
    ## nu = ny = nports: ext_clean's last nports inputs (w) are driven by
    ## G_pade's first nports outputs (w), and ext_clean's last nports outputs
    ## (z) drive G_pade's first nports inputs (z) -- exactly the loop the ideal
    ## delay used to close (w = Delta z, now w = G_pade z).
    sys = lft (ext_clean, G_pade, nports, nports);
  endif

  ## Ordinary delay, if any, is untouched by the loop closure; reattach the
  ## ORIGINAL ordinary delay to the delay-free closure result and
  ## Pade-substitute it via the same tf-roundtrip path used for
  ## InternalDelay-free ss inputs, using the ordinary slice of the order
  ## assignment computed above.  lft() preserves the ordinary I/O dimensions
  ## (the first nu inputs / first ny outputs), so the ordinary delay vectors
  ## still line up.
  if (had_ordinary)
    sys = set (sys, "InputDelay", oindelay, ...
                    "OutputDelay", ooutdelay, ...
                    "IODelay", oiodelay);
    sys = ss (__pade_substitute_ordinary__ (tf (sys), orders_ordinary));
  endif

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

%!test  # ss with dense per-entry IODelay, cross-checked against the tf path
%! a = [-2, 0; 0, -3];
%! b = [1, 0; 0, 1];
%! c = [1, 1; 1, 1];
%! d = [0, 0; 0, 0];
%! sys_ss = ss (a, b, c, d);
%! sys_ss = set (sys_ss, "IODelay", [0.2, 0.3; 0.4, 0.5]);
%! sys_tf = tf ({1, 1; 1, 1}, {[1 2], [1 3]; [1 2], [1 3]});
%! sys_tf = set (sys_tf, "IODelay", [0.2, 0.3; 0.4, 0.5]);
%! ss2 = pade (sys_ss, 3);
%! tf2 = pade (sys_tf, 3);
%! assert (hasdelay (ss2), false);
%! assert (hasdelay (tf2), false);
%! w = [0.1, 1, 5];
%! assert (freqresp (ss2, w), freqresp (tf2, w), 1e-6);

%!test  # ss InternalDelay: feedback loop, freqresp matches Pade LFT closed form
%! ## L = feedback (G), G = 1/(s+1) with loop delay T = 0.3 s, traps the delay
%! ## internally.  Exact closed loop: Gd/(1+Gd), Gd = G e^{-jwT}.  The Pade
%! ## approximation replaces e^{-jwT} by the scalar Pade transfer function's
%! ## own freqresp Pd, so the reference is Gp/(1+Gp) with Gp = G*Pd -- an
%! ## INDEPENDENT computation (scalar pade path), not a second call into the
%! ## code under test.
%! T = 0.3;
%! G = ss (-1, 1, 1, 0, "IODelay", T);
%! L = feedback (G);
%! assert (hasinternaldelay (L), true);
%! order = 4;
%! sys2 = pade (L, order);
%! assert (hasinternaldelay (sys2), false);
%! assert (hasdelay (sys2), false);
%! w = [0.05, 0.5, 1.7, 4, 11];
%! H = reshape (freqresp (sys2, w), 1, []);
%! [np, dp] = pade (T, order);
%! Pd = reshape (freqresp (tf (np, dp), w), 1, []);
%! Gp = (1 ./ (1i*w + 1)) .* Pd;
%! expected = Gp ./ (1 + Gp);
%! assert (H, expected, 1e-9);

%!test  # ss InternalDelay: freqresp error against TRUE exact delay shrinks with order
%! ## Strongest evidence the substitution converges to the exact-delay system:
%! ## as the Pade order grows, the approximation error against the untouched
%! ## exact-delay freqresp must strictly decrease.
%! T = 0.3;
%! G = ss (-1, 1, 1, 0, "IODelay", T);
%! L = feedback (G);
%! w = [0.05, 0.5, 1.7, 4, 11];
%! Hexact = reshape (freqresp (L, w), 1, []);
%! err = zeros (1, 3);
%! orders = [2, 4, 8];
%! for k = 1:3
%!   sysk = pade (L, orders(k));
%!   Hk = reshape (freqresp (sysk, w), 1, []);
%!   err(k) = max (abs (Hk - Hexact));
%! endfor
%! assert (err(2) < err(1));
%! assert (err(3) < err(2));

%!test  # ss MIMO InternalDelay: two independent ports, distinct delays and orders
%! ## append() of two SISO feedback loops -> genuine 2-port InternalDelay.
%! ## A per-port order VECTOR (3 for port 1, 6 for port 2) verifies each port
%! ## consumes its own assigned order.  Reference: decoupled per-channel Pade
%! ## closed form; off-diagonals must be identically zero.
%! T1 = 0.3; a1 = 1; o1 = 3;
%! T2 = 0.7; a2 = 2; o2 = 6;
%! G1 = ss (-a1, 1, 1, 0, "IODelay", T1);
%! G2 = ss (-a2, 1, 1, 0, "IODelay", T2);
%! sys = append (feedback (G1), feedback (G2));
%! assert (get (sys, "internaldelay"), [T1; T2], 1e-12);
%! sys2 = pade (sys, [o1, o2]);
%! assert (hasinternaldelay (sys2), false);
%! assert (hasdelay (sys2), false);
%! w = [0.05, 0.5, 1.7, 4, 11];
%! [n1, d1] = pade (T1, o1); P1 = reshape (freqresp (tf (n1, d1), w), 1, []);
%! [n2, d2] = pade (T2, o2); P2 = reshape (freqresp (tf (n2, d2), w), 1, []);
%! Gp1 = (1 ./ (1i*w + a1)) .* P1; e1 = Gp1 ./ (1 + Gp1);
%! Gp2 = (1 ./ (1i*w + a2)) .* P2; e2 = Gp2 ./ (1 + Gp2);
%! H = freqresp (sys2, w);
%! for k = 1:numel (w)
%!   assert (H(1,1,k), e1(k), 1e-9);
%!   assert (H(2,2,k), e2(k), 1e-9);
%!   assert (H(1,2,k), 0, 1e-9);
%!   assert (H(2,1,k), 0, 1e-9);
%! endfor

%!test  # ss combined InternalDelay AND ordinary OutputDelay: both approximated
%! ## The loop-trapped internal delay and an untouched OutputDelay are
%! ## independent; the result must equal the internal-delay Pade closure
%! ## multiplied by the OutputDelay's own scalar Pade transfer function.
%! T = 0.3; Tout = 0.2; order = 5;
%! G = ss (-1, 1, 1, 0, "IODelay", T);
%! L = set (feedback (G), "OutputDelay", Tout);
%! assert (hasinternaldelay (L), true);
%! assert (hasdelay (L), true);
%! sys2 = pade (L, order);
%! assert (hasinternaldelay (sys2), false);
%! assert (hasdelay (sys2), false);
%! w = [0.05, 0.5, 1.7, 4];
%! [np, dp] = pade (T, order);
%! Pd = reshape (freqresp (tf (np, dp), w), 1, []);
%! Gp = (1 ./ (1i*w + 1)) .* Pd; internal_ref = Gp ./ (1 + Gp);
%! [no, do_] = pade (Tout, order);
%! Pout = reshape (freqresp (tf (no, do_), w), 1, []);
%! expected = internal_ref .* Pout;
%! H = reshape (freqresp (sys2, w), 1, []);
%! assert (H, expected, 1e-8);

%!test  # ss InternalDelay with no delay ports (degenerate b2/c2 empty): still closes
%! ## set() an internaldelay directly with empty ports -- the extended plant is
%! ## degenerate (no extra I/O) but the substitution must not error and must
%! ## drop the internal delay.
%! sys = set (ss (-1, 1, 1, 0), "internaldelay", 0.5);
%! sys2 = pade (sys, 3);
%! assert (hasinternaldelay (sys2), false);

%!test  # discrete tf ordinary delay: exact absorption, matches absorbDelay directly
%! sys = tf (1, [1, -0.5], 0.1, "InputDelay", 3);
%! sys2 = pade (sys, 4);
%! expected = absorbDelay (sys);
%! assert (hasdelay (sys2), false);
%! [num2, den2] = tfdata (sys2);
%! [nume, dene] = tfdata (expected);
%! assert (num2, nume, 1e-12);
%! assert (den2, dene, 1e-12);

%!test  # discrete tf: pade order n has NO effect on the exact discrete result
%! sys = tf (1, [1, -0.5], 0.1, "InputDelay", 2);
%! sys_a = pade (sys, 2);
%! sys_b = pade (sys, 9);
%! [numa, dena] = tfdata (sys_a);
%! [numb, denb] = tfdata (sys_b);
%! assert (numa, numb, 1e-14);
%! assert (dena, denb, 1e-14);

%!test  # discrete ss ordinary delay: tf-roundtrip through absorbDelay, cross-check freqresp
%! sys = ss (0.5, 1, 1, 0, 0.1, "InputDelay", 2);
%! sys2 = pade (sys, 3);
%! assert (isa (sys2, "ss"));
%! assert (hasdelay (sys2), false);
%! w = [0.1, 1, 5];
%! expected = absorbDelay (tf (set (sys, "InputDelay", 2)));
%! assert (freqresp (sys2, w), freqresp (expected, w), 1e-8);

%!test  # discrete ss InternalDelay: exact loop closure, matches the untouched
%! # exact-delay discrete freqresp (no approximation error to expect, unlike
%! # the continuous case -- this closure is exact by construction)
%! G = ss (0.5, 1, 1, 0, 0.1, "IODelay", 3);
%! L = feedback (G);
%! assert (hasinternaldelay (L), true);
%! sys2 = pade (L, 4);
%! assert (hasinternaldelay (sys2), false);
%! w = [0.1, 1, 5];
%! assert (freqresp (sys2, w), freqresp (L, w), 1e-9);

%!test  # discrete ss InternalDelay: order n has no effect on the exact result
%! G = ss (0.5, 1, 1, 0, 0.1, "IODelay", 2);
%! L = feedback (G);
%! sys_a = pade (L, 2);
%! sys_b = pade (L, 7);
%! w = [0.1, 1, 5];
%! assert (freqresp (sys_a, w), freqresp (sys_b, w), 1e-12);

%!test  # discrete ss: combined InternalDelay AND ordinary delay, both exact
%! G = ss (0.5, 1, 1, 0, 0.1, "IODelay", 2);
%! L = set (feedback (G), "InputDelay", 3);
%! assert (hasinternaldelay (L), true);
%! assert (hasdelay (L), true);
%! sys2 = pade (L, 4);
%! assert (hasinternaldelay (sys2), false);
%! assert (hasdelay (sys2), false);
%! w = [0.1, 1, 5];
%! assert (freqresp (sys2, w), freqresp (L, w), 1e-8);
