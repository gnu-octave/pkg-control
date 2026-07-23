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
## @deftypefn {Function File} {@var{retsys} =} __sys_connect__ (@var{sys}, @var{M})
## This function is part of the Model Abstraction Layer.  No argument checking.
## For internal use only.
## @example
## @group
## Problem: Solve the system equations of
##   .
## E x(t) = A x(t) + B e(t)
##
##   y(t) = C x(t) + D e(t)
##
##   e(t) = u(t) + M y(t)
##
## in order to build
##   .
## K x(t) = F x(t) + G u(t)
##
##   y(t) = H x(t) + J u(t)
##
## Solution: Laplace Transformation
## E s X(s) = A X(s) + B U(s) + B M Y(s)                     [1]
##
##     Y(s) = C X(s) + D U(s) + D M Y(s)                     [2]
##
## solve [2] for Y(s)
## Y(s) = [I - D M]^(-1) C X(s)  +  [I - D M]^(-1) D U(s)
##
## substitute Z = [I - D M]^(-1)
## Y(s) = Z C X(s) + Z D U(s)                                [3]
##
## insert [3] in [1], solve for X(s)
## X(s) = [s E - (A + B M Z C)]^(-1) (B + B M Z D) U(s)      [4]
##
## inserting [4] in [3] finally yields
## Y(s) = Z C [s E - (A + B M Z C)]^(-1) (B + B M Z D) U(s)  +  Z D U(s)
##        \ /    |   \_____ _____/       \_____ _____/          \ /
##         H     K         F                   G                 J
## @end group
## @end example
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.3

function sys = __sys_connect__ (sys, m)

  ## FIXME: error for nonsensical stuff like  feedback (ss (1), ss (-1))
  ##        how to detect this?  all-zero descriptor matrix?  ! any (sys.e(:))
  ## TODO:  investigate whether all-zero sys.e make sense,
  ##        before disabling them in a rash decision.

  ## ----------------------------------------------------------------------
  ## Delay-aware dispatch.
  ##
  ## The grouped system handed to us may carry
  ##   (a) lti-level I/O delays  (InputDelay / OutputDelay / IODelay), and/or
  ##   (b) pre-existing internal delay ports  (b2/c2/d12/d21/d22, from an
  ##       earlier feedback/connect call carried through __sys_group__).
  ##
  ## When NEITHER is present -- the overwhelming common case -- we run the
  ## original, purely-rational loop-closure algebra below, unchanged, so that
  ## delay-free feedback()/connect() results are byte-for-byte identical to
  ## before this feature existed.
  ##
  ## Otherwise we hand off to __sys_connect_delay__ (see below), which first
  ## absorbs every nonzero I/O delay into an internal series delay port (an
  ## exact, transfer-preserving rewrite) and then closes the loop with a
  ## port-aware generalisation of the very same (I - D M)^{-1} algebra.
  ## ----------------------------------------------------------------------

  has_ports   = ! isempty (sys.b2);
  [id, od, iod] = get (sys, "inputdelay", "outputdelay", "iodelay");
  has_iodelay = any (id(:) != 0) || any (od(:) != 0) || any (iod(:) != 0);

  if (has_ports || has_iodelay)
    sys = __sys_connect_delay__ (sys, m, id, od, iod);
    return;
  endif

  a = sys.a;
  b = sys.b;
  c = sys.c;
  d = sys.d;

  z = eye (rows (d)) - d*m;

  if (rcond (z) >= eps)  # check for singularity
    
    sys.a = a + b*m/z*c;  # F
    sys.b = b + b*m/z*d;  # G
    sys.c = z\c;          # H
    sys.d = z\d;          # J
  
    ## sys.e remains constant: [] for ss models, e for dss models
    
  else
    
    ## construct descriptor model
    ## try to introduce the least states
    [pp, mm] = size (d);
    n = rows (a);
    if (isempty (sys.e))
      sys.e = eye (n);
    endif
    
    if (mm <= pp)
      ## Introduce state variable e = u + My
      ##   .
      ## E x =  A x +      B e + 0 u
      ##   .
      ## 0 e = MC x + (MD-I) e + I u  
      ##    
      ##   y =  C x +      D e + 0 u
      ##
      sys.a = [a, b; m*c, m*d-eye(mm)];
      sys.b = [zeros(n,mm); eye(mm)];
      sys.c = [c, d];
      sys.d = zeros (pp,mm);
      sys.e = blkdiag (sys.e, zeros (mm));
      sys.stname = [sys.stname; repmat({""}, mm, 1)];   
    else
      ## Introduce state variable y
      ##   .
      ## E x = A x + BM y + B u
      ##   .
      ## 0 y = C x -  Z y + D u  
      ##    
      ##   y = 0 x +  I y + 0 u
      ##          
      sys.a = [a, b*m; c, -z];
      sys.b = [b; d];
      sys.c = [zeros(pp, n), eye(pp)];
      sys.d = zeros (pp, mm);
      sys.e = blkdiag (sys.e, zeros (pp));
      sys.stname = [sys.stname; repmat({""}, pp, 1)];  
    endif

  endif

endfunction


## -*- texinfo -*-
## Delay-aware loop closure.  Called by __sys_connect__ only when the grouped
## system carries I/O delays or pre-existing internal delay ports.
##
## Strategy:
##
##  1. ABSORB every nonzero I/O delay of the grouped plant into an internal
##     series delay port.  A series delay port on input j is the exact rewrite
##
##         z = e_j ,   w = z(t - tau) ,   plant input j driven by w
##
##     (the direct e_j -> plant coupling b(:,j)/d(:,j) is moved onto the port),
##     and analogously on the output side.  This is a transfer-preserving
##     identity on the open-loop plant; it merely re-expresses an I/O time
##     shift as a (w = z(t-tau)) internal port so the subsequent purely
##     algebraic loop closure sees a delay-free rational plant.  The absorbed
##     lti delays are then zeroed.
##
##  2. CLOSE THE LOOP with a port-aware generalisation of the ordinary
##     Z = (I - D M)^{-1} algebra: the connection matrix M is now entirely
##     delay-free (all delays live in ports), so M is folded into the solve
##     exactly as before, while the delay ports (B2/C2/D12/D21/D22) are carried
##     through the same Z.
##
## SCOPE / partition decision:  only I/O delay on a channel that the
## connection matrix M actually references becomes an internal port; a
## delay on a channel M never touches (not part of any loop or cascade
## wiring) is left as an ordinary external I/O time-shift, unchanged.
## IODelay(i,j) is folded into input j's effective series delay, which is exact
## when input j feeds a single output (the block-diagonal-of-SISO topology that
## feedback()/connect() actually build).  Descriptor systems (nonempty E) and a
## singular delay-free connection (I - D M) with delays present are outside this
## scope and error clearly rather than return a wrong answer.
function sys = __sys_connect_delay__ (sys, m, id, od, iod)

  if (! isempty (sys.e))
    error ("__sys_connect__: internal delay with descriptor (E) models not supported");
  endif

  a = sys.a;
  b = sys.b;
  c = sys.c;
  d = sys.d;

  n = rows (a);
  [p, mm] = size (d);

  ## ---- initialise delay-port blocks (existing ports, or empty) ----------
  if (isempty (sys.b2))
    B2  = zeros (n, 0);
    C2  = zeros (0, n);
    D12 = zeros (p, 0);
    D21 = zeros (0, mm);
    D22 = zeros (0, 0);
    tau = zeros (0, 1);
  else
    B2  = sys.b2;
    C2  = sys.c2;
    D12 = sys.d12;
    D21 = sys.d21;
    D22 = sys.d22;
    tau = sys.internaldelay(:);
  endif

  ## ---- 1. absorb I/O delays into series ports ---------------------------
  ## effective per-input series delay: InputDelay(j) + column j of IODelay
  ##
  ## v(i,j) = InputDelay(j) + IODelay(i,j) is the delay row i would need
  ## input j's series port to carry.  Only rows with a genuinely REQUESTED
  ## (nonzero) IODelay entry are compared for agreement -- rows with
  ## IODelay(i,j) == 0 make no request at all (this also naturally excludes
  ## rows belonging to an unrelated block padded in by __sys_group__, e.g.
  ## the identity compensator feedback() appends, whose IODelay is always
  ## zero).  Folding column j into a single shared port is exact exactly
  ## when every REQUESTING row agrees on the value: the common case (plain
  ## InputDelay, or a diagonal IODelay column) and the equal-nonzero-value
  ## case (e.g. IODelay = [0.3; 0.3]), which used to be mis-summed to 0.6
  ## and then rejected by the old any(sum(iod!=0,1)>1) guard below.  When
  ## requesting rows genuinely disagree in value, a single shared
  ## state-space column b(:,j) cannot carry two different delays to two
  ## different outputs; that decomposition is not implemented yet (error).
  inport = id(:).';                       # 1-by-m
  for j = 1 : mm
    vals = unique (iod (iod(:, j) != 0, j));
    if (numel (vals) > 1)
      error (["__sys_connect__: an input feeding multiple differently-delayed ", ...
              "outputs through IODelay is not yet supported"]);
    elseif (numel (vals) == 1)
      inport(j) += vals;
    endif
  endfor
  outport = od(:).';                      # 1-by-p

  ## A channel's I/O delay is promoted into an internal series delay port ONLY
  ## when that channel actually participates in the connection matrix M (i.e.
  ## it is genuinely cascaded or fed back).  A delay on a channel that M never
  ## touches -- e.g. the free external InputDelay of sys1 in a feedforward
  ## sys2*sys1 cascade, whose input row of M is all zero -- is NOT a trapped
  ## internal delay: it must pass straight through as an ordinary I/O time
  ## shift, exactly as before internal-delay support existed.  Absorbing it
  ## unconditionally would force a spurious internal-delay realisation (and
  ## block tf/zpk conversion) for a plain cascade.  M has one row per input
  ## channel and one column per output channel (e = u + M y).
  in_touched  = any (m != 0, 2).';        # 1-by-m: input row j has a nonzero
  out_touched = any (m != 0, 1);          # 1-by-p: output column i has a nonzero

  absorbed_in  = false (1, mm);
  absorbed_out = false (1, p);

  ## input-side ports first (they append columns that output ports may read)
  for j = 1 : mm
    if (inport(j) > 0 && in_touched(j))
      q = columns (B2);
      B2  = [B2,  b(:,j)];
      D12 = [D12, d(:,j)];
      C2  = [C2;  zeros(1,n)];
      D21 = [D21; full(sparse(1,j,1,1,mm))];
      D22 = [D22, zeros(q,1); zeros(1,q), 0];
      tau = [tau; inport(j)];
      b(:,j) = 0;                         # input j now reaches plant via port
      d(:,j) = 0;
      absorbed_in(j) = true;
    endif
  endfor

  ## output-side ports (read the possibly-updated c/d/D12 rows)
  for i = 1 : p
    if (outport(i) > 0 && out_touched(i))
      q = columns (B2);
      ## new port picks up the internal (pre-delay) output-i signal:
      ##   z = c(i,:) x + d(i,:) e + D12(i,:) w_existing
      B2  = [B2,  zeros(n,1)];            # port does not drive states
      C2  = [C2;  c(i,:)];
      D21 = [D21; d(i,:)];
      D22 = [D22, zeros(q,1); D12(i,:), 0];
      ## external output i is now emitted through the delayed port only
      c(i,:)   = 0;
      d(i,:)   = 0;
      D12(i,:) = 0;
      newcol   = zeros (p, 1);
      newcol(i) = 1;                      # external output i = new (delayed) w
      D12 = [D12, newcol];
      tau = [tau; outport(i)];
      absorbed_out(i) = true;
    endif
  endfor

  ## Clear ONLY the I/O delays that were actually absorbed into ports; delays
  ## on channels that M never touched pass through unchanged.  IODelay(i,j) is
  ## folded into input j's series delay, so it is cleared exactly when input j
  ## was absorbed.
  id_new  = id(:);   id_new(absorbed_in)   = 0;
  od_new  = od(:);   od_new(absorbed_out)  = 0;
  iod_new = iod;     iod_new(:, absorbed_in) = 0;
  sys = set (sys, "inputdelay", id_new, ...
                  "outputdelay", od_new, ...
                  "iodelay", iod_new);

  ## ---- 2. port-aware algebraic loop closure  e = u + M y ----------------
  z = eye (p) - d*m;

  if (rcond (z) < eps)
    error ("__sys_connect__: singular algebraic loop with internal delays not supported");
  endif

  sys.a   = a + b*m/z*c;                  # F
  sys.b   = b + b*m/z*d;                  # G
  sys.c   = z\c;                          # H
  sys.d   = z\d;                          # J

  sys.b2  = B2  + b*m/z*D12;
  sys.c2  = C2  + D21*m/z*c;
  sys.d12 = z\D12;
  sys.d21 = D21 + D21*m/z*d;
  sys.d22 = D22 + D21*m/z*D12;
  sys.internaldelay = tau;

  sys.e = [];

endfunction


## Regression: a delay on a channel the connection matrix M never touches must
## pass straight through as an ordinary I/O delay -- NOT be promoted into an
## internal delay port.  In a feedforward cascade sys2*sys1, sys1's InputDelay
## is a free external input (its input row of M is all zero), so the mtimes
## result must keep it as InputDelay and carry no internal delay.  This is the
## exact situation that made @lti/series.m's MIMO-guard test spuriously fail
## with a sys2tf/InternalDelay error before this narrowing was applied.
%!test
%! s1 = ss (-1, 1, 1, 0, "InputDelay", 0.3);
%! s2 = ss (-2, 1, 1, 0);
%! sys = s2 * s1;                          # feedforward cascade -> mtimes
%! assert (get (sys, "inputdelay"), 0.3, 1e-12);
%! assert (isempty (get (sys, "internaldelay")));


## Bug fix: a column with 2+ EQUAL nonzero IODelay entries must absorb into
## a single shared port of that value, not error and not silently mis-sum
## (id(j) + sum(iod(:,j)) used to double-count equal values, e.g. 0.3+0.3 =
## 0.6, which is why the old any(sum(iod!=0,1)>1) guard rejected this
## topology outright rather than risk returning that wrong answer).
%!test
%! a = -eye (2);
%! b = eye (2);
%! c = eye (2);
%! d = zeros (2);
%! iod = [0.3, 0; 0.3, 0];                 # column 1: two EQUAL nonzero entries
%! G = ss (a, b, c, d, "IODelay", iod);
%! cl = feedback (G);                      # unity negative feedback, both channels
%! assert (get (cl, "internaldelay"), 0.3, 1e-10);   # single 0.3, not 0.6
%! ## hand-derived closed-form: the plant is decoupled (a diagonal, b = c =
%! ## I), so each channel is an independent SISO loop; channel 1 carries the
%! ## series delay e^{-0.3s}, channel 2 has none.
%! ##   T1(s) = g(s) e^{-0.3s} / (1 + g(s) e^{-0.3s}),   g(s) = 1/(s+1)
%! ##   T2(s) = g(s) / (1 + g(s))
%! w = 1;
%! g = 1 / (1i*w + 1);
%! T1 = g * exp (-1i*0.3*w) / (1 + g * exp (-1i*0.3*w));
%! T2 = g / (1 + g);
%! fr = freqresp (cl, w);
%! assert (fr(:, :, 1), [T1, 0; 0, T2], 1e-8);


## Exotic topology regression: nested feedback around a delayed inner loop.
## The inner loop's IODelay is trapped by the inner feedback() call, then the
## inner closed loop is cascaded with a further controller and closed again
## by an outer feedback() call.  This already worked before this task (the
## inner delay is a single-input/single-output series port throughout, so
## the uniform-column path always applied); promoted here into a permanent
## regression test with a hand-derived closed-form reference.
%!test
%! P = ss (-1, 1, 1, 0, "IODelay", 0.2);
%! H = ss (-3, 1, 1, 0);
%! Inner = feedback (P, H);                # inner loop traps the delay
%! assert (get (Inner, "internaldelay"), 0.2, 1e-10);
%! Cc = ss (-2, 1, 1, 0);
%! Outer = feedback (Inner * Cc);          # outer unity feedback
%! assert (get (Outer, "internaldelay"), 0.2, 1e-10);
%! ## hand derivation:  p(s) = 1/(s+1), h(s) = 1/(s+3), c(s) = 1/(s+2)
%! ##   Inner(s) = p e^{-0.2s} / (1 + p e^{-0.2s} h)
%! ##   Outer(s) = Inner*c / (1 + Inner*c)
%! w = 1.3;
%! p = 1 / (1i*w + 1);
%! h = 1 / (1i*w + 3);
%! inner_s = p * exp (-1i*0.2*w) / (1 + p * exp (-1i*0.2*w) * h);
%! c = 1 / (1i*w + 2);
%! outer_s = inner_s * c / (1 + inner_s * c);
%! assert (freqresp (Outer, w), outer_s, 1e-8);


## Exotic topology regression: a 3-block delay ring closed via connect(),
## each block carrying its own distinct IODelay.  Already worked before this
## task (each block's delay is absorbed as its own single-input series
## port); promoted here into a permanent regression test.
%!test
%! G1 = ss (-1, 1, 1, 0, "IODelay", 0.1);
%! G2 = ss (-2, 1, 1, 0, "IODelay", 0.15);
%! G3 = ss (-3, 1, 1, 0, "IODelay", 0.05);
%! Gall = append (G1, G2, G3);
%! ## ring:  u1 = r - y3 ,  u2 = y1 ,  u3 = y2
%! cm = [1, -3; 2, 1; 3, 2];
%! Ring = connect (Gall, cm, 1, 1);
%! assert (sort (get (Ring, "internaldelay")), sort ([0.1; 0.15; 0.05]), 1e-10);
%! ## hand derivation:  g_k(s) = 1/(s+k) e^{-tau_k s},  loop gain L = g1 g2 g3
%! ##   Ring(s) = g1 / (1 + L)         (transfer from external input r to y1)
%! w = 0.7;
%! g1 = 1 / (1i*w + 1) * exp (-1i*0.1*w);
%! g2 = 1 / (1i*w + 2) * exp (-1i*0.15*w);
%! g3 = 1 / (1i*w + 3) * exp (-1i*0.05*w);
%! ring_s = g1 / (1 + g1*g2*g3);
%! assert (freqresp (Ring, w), ring_s, 1e-8);


## Exotic topology regression: cross-coupled two-block feedback where each
## block carries its own distinct IODelay (as opposed to a single scalar
## delay in only the forward path).  Already worked before this task;
## promoted here into a permanent regression test.
%!test
%! Ga = ss (-1, 1, 1, 0, "IODelay", 0.2);
%! Gb = ss (-2, 1, 1, 0, "IODelay", 0.4);
%! Cross = feedback (Ga, Gb);
%! assert (sort (get (Cross, "internaldelay")), sort ([0.2; 0.4]), 1e-10);
%! ## hand derivation:  ga(s) = 1/(s+1) e^{-0.2s},  gb(s) = 1/(s+2) e^{-0.4s}
%! ##   Cross(s) = ga / (1 + ga*gb)
%! w = 0.9;
%! ga = 1 / (1i*w + 1) * exp (-1i*0.2*w);
%! gb = 1 / (1i*w + 2) * exp (-1i*0.4*w);
%! cross_s = ga / (1 + ga*gb);
%! assert (freqresp (Cross, w), cross_s, 1e-8);
