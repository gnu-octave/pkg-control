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
## @deftypefn {Function File} {@var{sys} =} c2d (@var{sys}, @var{tsam})
## @deftypefnx {Function File} {@var{sys} =} c2d (@var{sys}, @var{tsam}, @var{method})
## @deftypefnx {Function File} {@var{sys} =} c2d (@var{sys}, @var{tsam}, @var{'prewarp'}, @var{w0})
## Convert the continuous @acronym{LTI} model into its discrete-time equivalent.
##
## @strong{Inputs}
## @table @var
## @item sys
## Continuous-time @acronym{LTI} model.
## @item tsam
## Sampling time in seconds.
## @item method
## Optional conversion method.  If not specified, default method @var{"zoh"}
## is taken.
## @table @var
## @item 'impulse'
## Impulse Invarient transformation.
## @item 'zoh'
## Zero-order hold or matrix exponential.
## @item 'foh'
## First-order hold, linear approximation of the input signals between
## two sample times
## @item 'tustin', 'bilin'
## Bilinear transformation or Tustin approximation.
## @item 'prewarp'
## Bilinear transformation with pre-warping at frequency @var{w0}.
## @item 'matched'
## Matched pole/zero method.
## @end table
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## Discrete-time @acronym{LTI} model.
## @end table
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.3

function sys = c2d (sys, tsam, method = "std", w0 = 0)

  if (nargin < 2 || nargin > 4)
    print_usage ();
  endif

  if (! isa (sys, "lti"))
    error ("c2d: first argument is not an LTI model");
  endif

  if (isdt (sys))
    error ("c2d: system is already discrete-time");
  endif

  if (! issample (tsam))
    error ("c2d: second argument is not a valid sample time");
  endif

  delay_modeling = "delay";

  if (isstruct (method))
    if (w0 != 0)
      error ("c2d: cannot combine a c2dOptions struct with an explicit prewarp frequency argument");
    endif

    opt = method;
    method = opt.Method;
    w0 = opt.PrewarpFrequency;
    delay_modeling = opt.DelayModeling;
  endif

  if (! ischar (method))
    error ("c2d: third argument is not a string");
  endif

  if (! issample (w0, 0))
    error ("c2d: fourth argument is not a valid pre-warping frequency");
  endif

  origsys = sys;

  if (hasinternaldelay (origsys))
    ## Discretizing an InternalDelay system is exactly discretizing the
    ## extended ordinary system (A, [B1 B2], [C1;C2], [[D11 D12];[D21 D22]])
    ## via the unchanged __c2d__ (which only looks at matrix shapes), then
    ## rounding tau to whole samples the same way InputDelay/OutputDelay are
    ## rounded below.
    ## __ss_ext_split__ copies InputDelay/OutputDelay/IODelay unchanged from
    ## origsys (still in seconds), so if origsys also carries ordinary I/O
    ## delay, round it to samples here, the same way the hasdelay() branch
    ## below does for the InternalDelay-free case.  Combining InternalDelay
    ## with DelayModeling="state" or fractional (Thiran) ordinary delay is
    ## out of scope; only the plain round-to-samples case is handled.
    [ext_sys, nu, ny] = __ss_ext_build__ (origsys);
    ext_sys = __c2d__ (ext_sys, tsam, lower (method), w0);
    sys = __ss_ext_split__ (origsys, ext_sys, nu, ny);
    sys.tsam = tsam;
    tau = get (origsys, "internaldelay");
    tau_samples = round (tau / tsam);

    ## A nonzero InternalDelay that rounds to 0 samples would silently drop
    ## that port's feedthrough during simulation instead of solving the
    ## resulting algebraic loop (documented in __delay_lookup__/
    ## __buffered_sim__'s tau>=1 assumption) -- convert this from a
    ## silent-wrong-answer risk into a clear error, consistent with every
    ## other guard on this branch.
    if (any (tau != 0 & tau_samples == 0))
      bad = tau (find (tau != 0 & tau_samples == 0, 1));
      error (["c2d: InternalDelay of %g seconds is too small relative to ", ...
              "the sampling time %g and would round to 0 samples"], ...
             bad, tsam);
    endif

    sys = set (sys, "internaldelay", tau_samples);

    if (hasdelay (origsys))
      [indelay, outdelay, iodelay] = get (origsys, "inputdelay", "outputdelay", "iodelay");
      sys = set (sys, "InputDelay", round (indelay / tsam), ...
                      "OutputDelay", round (outdelay / tsam), ...
                      "IODelay", round (iodelay / tsam));
    endif

    return;
  endif

  if (hasdelay (origsys))
    [indelay, outdelay, iodelay] = get (origsys, "inputdelay", "outputdelay", "iodelay");

    indelay_samples = indelay / tsam;
    outdelay_samples = outdelay / tsam;
    iodelay_samples = iodelay / tsam;

    all_samples = [indelay_samples(:); outdelay_samples(:); iodelay_samples(:)];

    if (strcmp (delay_modeling, "state"))
      % Delay must be zeroed before __c2d__ for the same recursive-re-entry
      % reason documented below (ss() preserves delay fields), even though
      % this branch never calls thiran -- otherwise a MIMO/recursive call
      % into __c2d__ would re-enter this same "state" branch again.
      sys_no_delay = set (origsys, "InputDelay", 0, "OutputDelay", 0, "IODelay", 0);

      % absorbDelay errors on ss -- route ss inputs through tf, absorb
      % there, then convert back, since absorbDelay's per-channel math
      % already handles a general IODelay matrix correctly.
      is_ss_input = isa (origsys, "ss");

      if (is_ss_input)
        sys = tf (sys_no_delay);
      else
        sys = sys_no_delay;
      endif

      sys = __c2d__ (sys, tsam, lower (method), w0);
      sys.tsam = tsam;
      sys = set (sys, "InputDelay", round (indelay_samples), ...
                      "OutputDelay", round (outdelay_samples), ...
                      "IODelay", round (iodelay_samples));
      sys = absorbDelay (sys);

      if (is_ss_input)
        sys = ss (sys);
      endif
    elseif (any (abs (all_samples - round (all_samples)) > 1e-8))
      [p, m] = size (origsys);

      if (any (strcmpi (method, {"zoh", "std"})) && __c2d_frac_applicable__ (origsys, tsam))
        sys = __c2d_frac_ss__ (origsys, tsam);
        sys.tsam = tsam;
      elseif (p != 1 || m != 1)
        error ("c2d: fractional delays on MIMO systems are not yet supported (per-channel Thiran approximation is not yet implemented)");
      else
        % __c2d__ recurses into this same function via ss() for the general
        % case, and ss() does not strip delay fields -- so discretizing
        % origsys directly here would re-enter this fractional branch a
        % second time and double-apply the Thiran filter below.  Delay must
        % be zeroed before the __c2d__ call, not just after.
        sys_no_delay = set (origsys, "InputDelay", 0, "OutputDelay", 0, "IODelay", 0);
        sys = __c2d__ (sys_no_delay, tsam, lower (method), w0);
        sys.tsam = tsam;

        total_delay = totaldelay (origsys);
        filt = thiran (total_delay, tsam);
        sys = sys * filt;
        sys = set (sys, "InputDelay", 0, "OutputDelay", 0, "IODelay", 0);
      endif
    else
      sys = __c2d__ (origsys, tsam, lower (method), w0);
      sys.tsam = tsam;
      sys = set (sys, "InputDelay", round (indelay_samples), ...
                      "OutputDelay", round (outdelay_samples), ...
                      "IODelay", round (iodelay_samples));
    endif
  else
    sys = __c2d__ (origsys, tsam, lower (method), w0);
    sys.tsam = tsam;
  endif

endfunction


## bilinear transformation
## using oct-file directly
%!shared Mo, Me
%! A = [  1.0  0.5
%!        0.5  1.0 ].';
%!
%! B = [  0.0 -1.0
%!        1.0  0.0 ].';
%!
%! C = [ -1.0  0.0
%!        0.0  1.0 ].';
%!
%! D = [  1.0  0.0
%!        0.0 -1.0 ].';
%!
%! [Ao, Bo, Co, Do] = __sl_ab04md__ (A, B, C, D, 1.0, 1.0, false);
%!
%! Ae = [ -1.0000  -4.0000
%!        -4.0000  -1.0000 ];
%!
%! Be = [  2.8284   0.0000
%!         0.0000  -2.8284 ];
%!
%! Ce = [  0.0000   2.8284
%!        -2.8284   0.0000 ];
%!
%! De = [ -1.0000   0.0000
%!         0.0000  -3.0000 ];
%!
%! Mo = [Ao, Bo; Co, Do];
%! Me = [Ae, Be; Ce, De];
%!
%!assert (Mo, Me, 1e-4);


## bilinear transformation
## user function
%!shared Mo, Me
%! A = [  1.0  0.5
%!        0.5  1.0 ].';
%!
%! B = [  0.0 -1.0
%!        1.0  0.0 ].';
%!
%! C = [ -1.0  0.0
%!        0.0  1.0 ].';
%!
%! D = [  1.0  0.0
%!        0.0 -1.0 ].';
%!
%! [Ao, Bo, Co, Do] = ssdata (c2d (ss (A, B, C, D), 2, "tustin"));
%!
%! Ae = [ -1.0000  -4.0000
%!        -4.0000  -1.0000 ];
%!
%! Be = [  2.8284   0.0000
%!         0.0000  -2.8284 ];
%!
%! Ce = [  0.0000   2.8284
%!        -2.8284   0.0000 ];
%!
%! De = [ -1.0000   0.0000
%!         0.0000  -3.0000 ];
%!
%! Mo = [Ao, Bo; Co, Do];
%! Me = [Ae, Be; Ce, De];
%!
%!assert (Mo, Me, 1e-4);


## impulse invariant
%!shared Mo, Me
%! A = [  1.0  0.5
%!        0.5  1.0 ];
%!
%! B = [  0.0 -1.0
%!        1.0  0.0 ];
%!
%! C = [ -1.0  0.0
%!        0.0  1.0 ];
%!
%! D = [  0.0  0.0
%!        0.0  0.0 ];
%!
%! [Ao, Bo, Co, Do] = ssdata (c2d (ss(A,B,C,D), 2, "imp"));
%!
%! Ae = [ 11.4019   8.6836
%!         8.6836  11.4019 ];
%!
%! Be = [ 17.3673 -22.8038
%!        22.8038 -17.3673 ];
%!
%! Ce = [ -1.0000   0.0000
%!         0.0000   1.0000 ];
%!
%! De = [  0.0000   2.0000
%!         2.0000   0.0000 ];
%!
%! Mo = [Ao, Bo; Co, Do];
%! Me = [Ae, Be; Ce, De];
%!
%!assert (Mo, Me, 1e-4);


## impulse invariant for transfer function
%!shared Mo, Me
%! G = tf ({[1 0],1;[1],1},{[1 1 1],[1 1];[1 0],[1 2 1]});
%!
%! [nuo, dno] = tfdata (c2d (G, 2, "imp"));
%!
%! nue = {[2 -0.3011 0], [2 0]; [2 0], [0 0.5413 0]};
%! dne = {[1 0.1181 0.1353], [1 -0.1353]; [1 -1], [1 -0.2707 0.01832]};
%!
%! Mo = [ nuo{1,1} nuo{1,2} nuo{2,1} nuo{2,2} dno{1,1} dno{1,2} dno{2,1} dno{2,2} ];
%! Me = [ nue{1,1} nue{1,2} nue{2,1} nue{2,2} dne{1,1} dne{1,2} dne{2,1} dne{2,2} ];
%!
%!assert (Mo, Me, 1e-4);


## bilinear transformation
## both directions
%!shared Mo, Me
%! A = [  1.0  0.5
%!        0.5  1.0 ];
%!
%! B = [  0.0 -1.0
%!        1.0  0.0 ];
%!
%! C = [ -1.0  0.0
%!        0.0  1.0 ];
%!
%! D = [  1.0  0.0
%!        0.0 -1.0 ];
%!
%! [Ao, Bo, Co, Do] = ssdata (d2c (c2d (ss (A, B, C, D), 2, "tustin"), "tustin"));
%!
%! Mo = [Ao, Bo; Co, Do];
%! Me = [A, B; C, D];
%!
%!assert (Mo, Me, 1e-4);


## zero-order hold
## both directions
%!shared Mo, Me
%! A = [  1.0  0.5
%!        0.5  1.0 ];
%!
%! B = [  0.0 -1.0
%!        1.0  0.0 ];
%!
%! C = [ -1.0  0.0
%!        0.0  1.0 ];
%!
%! D = [  1.0  0.0
%!        0.0 -1.0 ];
%!
%! [Ao, Bo, Co, Do] = ssdata (d2c (c2d (ss (A, B, C, D), 2, "zoh"), "zoh"));
%!
%! Mo = [Ao, Bo; Co, Do];
%! Me = [A, B; C, D];
%!
%!assert (Mo, Me, 1e-4);


## first-order hold
## user function
%!shared Mo, Me
%! A = [  1.0  0.5
%!        0.5  1.0 ];
%!
%! B = [  0.0 -1.0
%!        1.0  0.0 ];
%!
%! C = [ -1.0  0.0
%!        0.0  1.0 ];
%!
%! D = [  1.0  0.0
%!        0.0 -1.0 ];
%!
%! [Ao, Bo, Co, Do] = ssdata (c2d (ss (A, B, C, D), 2, "foh"));
%!
%! Ae = [ 11.4019    8.6836
%!         8.6836   11.4019 ];
%!
%! Be = [  37.5206  -43.4256
%!         43.4256  -37.5206 ];
%!
%! Ce = [ -1.0000   0.0000
%!         0.0000   1.0000 ];
%!
%! De = [ -0.0690   2.5056
%!         2.5056  -2.0690 ];
%!
%! Mo = [Ao, Bo; Co, Do];
%! Me = [Ae, Be; Ce, De];
%!
%!assert (Mo, Me, 1e-4);


## bilinear transformation with pre-warping
## both directions
%!shared Mo, Me
%! A = [  1.0  0.5
%!        0.5  1.0 ];
%!
%! B = [  0.0 -1.0
%!        1.0  0.0 ];
%!
%! C = [ -1.0  0.0
%!        0.0  1.0 ];
%!
%! D = [  1.0  0.0
%!        0.0 -1.0 ];
%!
%! [Ao, Bo, Co, Do] = ssdata (d2c (c2d (ss (A, B, C, D), 2, "prewarp", 1000), "prewarp", 1000));
%!
%! Mo = [Ao, Bo; Co, Do];
%! Me = [A, B; C, D];
%!
%!assert (Mo, Me, 1e-4);


## matrix exponential
%!shared Aex, Aexint, Aex_exp, Aexint_exp
%! A = [  5.0   4.0   3.0   2.0   1.0
%!        1.0   6.0   0.0   4.0   3.0
%!        2.0   0.0   7.0   6.0   5.0
%!        1.0   3.0   1.0   8.0   7.0
%!        2.0   5.0   7.0   1.0   9.0 ];
%!
%! Aex_exp = [ 1.8391   0.9476   0.7920   0.8216   0.7811
%!             0.3359   2.2262   0.4013   1.0078   1.0957
%!             0.6335   0.6776   2.6933   1.6155   1.8502
%!             0.4804   1.1561   0.9110   2.7461   2.0854
%!             0.7105   1.4244   1.8835   1.0966   3.4134 ];
%!
%! Aexint_exp = [ 0.1347   0.0352   0.0284   0.0272   0.0231
%!                0.0114   0.1477   0.0104   0.0369   0.0368
%!                0.0218   0.0178   0.1624   0.0580   0.0619
%!                0.0152   0.0385   0.0267   0.1660   0.0732
%!                0.0240   0.0503   0.0679   0.0317   0.1863 ];
%!
%! [Aex, Aexint] = __sl_mb05nd__ (A, 0.1, 0.0001);
%!
%!assert (Aex, Aex_exp, 1e-4);
%!assert (Aexint, Aexint_exp, 1e-4);


%!test  # exact integer-sample InputDelay survives c2d
%! sys = tf (1, [1 1], "InputDelay", 1.0);
%! dsys = c2d (sys, 0.5, "zoh");
%! assert (isdt (dsys), true);
%! assert (dsys.InputDelay, 2);
%! assert (dsys.OutputDelay, 0);
%! assert (dsys.IODelay, 0);

%!test  # no delay: unaffected (regression)
%! sys = tf (1, [1 1]);
%! dsys = c2d (sys, 0.5, "zoh");
%! assert (hasdelay (dsys), false);


%!test  # SISO tf fractional InputDelay: exact split-ZOH via zoh (default)
%! sys = tf (1, [1 1], "InputDelay", 1.33);
%! dsys = c2d (sys, 0.5, "zoh");
%! assert (isdt (dsys), true);
%! ## d = floor(tau/tsam) = floor(1.33/0.5) = 2, plus 1 more because the
%! ## fractional remainder (0.33 s) is nonzero -- same convention as
%! ## __c2d_frac_zoh__/__c2d_frac_ss__.
%! assert (hasdelay (dsys), true);
%! assert (dsys.InputDelay, 3);
%! [numd, dend] = tfdata (dsys, "v");
%! ## independent ground truth: split-ZOH construction (Astrom & Wittenmark),
%! ## computed directly via expm/integral (not via the code under test) --
%! ## see the reproduction script in the task report for how these numbers
%! ## were derived: A=-1, B=1, tsam=0.5, tau=1.33, d=floor(tau/tsam)=2,
%! ## tau_frac=0.33.
%! nume = [0.1563351834, 0.2371341569];
%! dene = [1, -0.6065306597];
%! assert (numd, nume, 1e-6);
%! assert (dend, dene, 1e-6);
%!
%! ## Thiran approximation path is still reachable for a non-zoh method.
%! dsys_tustin = c2d (sys, 0.5, "tustin");
%! assert (hasdelay (dsys_tustin), false);

%!test  # SISO zpk fractional OutputDelay: exact split-ZOH via zoh (default)
%! sys = zpk ([], -1, 1, "OutputDelay", 1.33);
%! dsys = c2d (sys, 0.5, "zoh");
%! assert (isdt (dsys), true);
%! assert (hasdelay (dsys), true);
%! assert (dsys.InputDelay, 3);
%! w = [0.1, 1, 5];
%! ## independent ground truth: same split-ZOH construction as the tf
%! ## InputDelay test above -- the dynamics (pole at -1, gain 1) are
%! ## identical, since __c2d_frac_ss__ collapses all delay fields into a
%! ## scalar totaldelay () before splitting, so the same reference numbers
%! ## apply.
%! nume = [0.1563351834, 0.2371341569];
%! dene = [1, -0.6065306597];
%! tsam = 0.5;
%! expected = tf (nume, dene, tsam);
%! expected = set (expected, "InputDelay", 3);
%! assert (freqresp (dsys, w), freqresp (expected, w), 1e-6);
%!
%! ## Thiran approximation path is still reachable for a non-zoh method.
%! dsys_tustin = c2d (sys, 0.5, "tustin");
%! assert (hasdelay (dsys_tustin), false);

%!test  # SIMO (single input, multiple outputs) fractional InputDelay:
%! ## Finding 1 fix-pass regression test.  Before the fix, __c2d_frac_ss__'s
%! ## trailing tf/zpk pole/zero-at-origin cancellation block was guarded
%! ## only by "m == 1", which is trivially true here (m == 1) even though
%! ## the system is not truly SISO (p == 2) -- the cancellation code calls
%! ## tfdata(sys, "v") in a way that assumes a genuinely SISO system, and
%! ## crashed with "ERROR: abs: not defined for cell".  The guard is now
%! ## "p == 1 && m == 1", so this case is correctly routed through the
%! ## general per-column MIMO augmentation path instead (no cancellation
%! ## attempted).  Verify no crash, and cross-check each output channel's
%! ## response against an independently, separately discretized SISO
%! ## channel.
%! sys = tf ({1;1}, {[1 1];[1 2]}, "InputDelay", 1.33);
%! dsys = c2d (sys, 0.5, "zoh");
%! assert (isdt (dsys), true);
%! h1 = tf (1, [1 1], "InputDelay", 1.33);
%! h2 = tf (1, [1 2], "InputDelay", 1.33);
%! d1 = c2d (h1, 0.5, "zoh");
%! d2 = c2d (h2, 0.5, "zoh");
%! w = [0.1, 1, 5];
%! resp = freqresp (dsys, w);
%! assert (reshape (resp(1,1,:), 1, []), reshape (freqresp (d1, w), 1, []), 1e-10);
%! assert (reshape (resp(2,1,:), 1, []), reshape (freqresp (d2, w), 1, []), 1e-10);

%!test  # MIMO fractional InputDelay, uniform per column: now handled exactly
%! ## by the Task 2 per-column split-ZOH generalization (this used to error
%! ## before MIMO support was added -- __c2d_frac_applicable__ now accepts
%! ## it because each column's total delay is constant down that column).
%! ## InputDelay is [2;0], not [3;0]: the fix-pass correction removed the
%! ## blanket "+1" applied to every fractional-delay column, since the
%! ## extra register state remains part of the augmented dynamics here (no
%! ## tf/zpk pole-zero cancellation happens for p>1) and its own
%! ## eigenvalue-zero pole already contributes exactly the missing sample.
%! sys = tf ({1,1;1,1}, {[1 1],[1 2];[1 3],[1 4]}, "InputDelay", [1.33;0]);
%! dsys = c2d (sys, 0.5, "zoh");
%! assert (isdt (dsys), true);
%! assert (dsys.InputDelay, [2;0]);

%!test  # MIMO: column 1 exact-integer (nonzero) delay, column 2 fractional --
%! ## only column 2 gets an extra register state; column 1's response must
%! ## match a plain (no-extra-state) single-channel zoh discretization
%! ## exactly, confirming the exact-integer column is untouched by the
%! ## per-column augmentation applied to its fractional sibling.
%! A = [-1, 0; 0, -3];
%! B = [1, 0; 0, 1];
%! C = [1, 0; 0, 1];
%! D = [0, 0; 0, 0];
%! sys = ss (A, B, C, D);
%! sys = set (sys, "InputDelay", [0.2; 0.37]);   # col1: exact 2 samples @ Ts=0.1
%! tsam = 0.1;
%! dsys = c2d (sys, tsam);
%! assert (dsys.InputDelay(1), 2);
%! assert (dsys.InputDelay(2), 3);   # floor(0.37/0.1)=3; no +1, since the
%!                                    # register survives as a real state in
%!                                    # this ss/MIMO output (only cancelled-out
%!                                    # tf/zpk SISO conversions get the +1)
%! sys1 = ss (A(1,1), B(1,1), C(1,1), D(1,1));
%! sys1 = set (sys1, "InputDelay", 0.2);
%! dsys1 = c2d (sys1, tsam);
%! w = [0.1, 1, 5];
%! resp_full = freqresp (dsys, w);
%! resp_ref = freqresp (dsys1, w);
%! assert (squeeze (resp_full(1,1,:)), squeeze (resp_ref), 1e-8);

%!test  # MIMO fractional OutputDelay, uniform across all outputs: in-scope
%! ## (Finding 2 fix-pass): OutputDelay may be nonzero for MIMO as long as
%! ## it is uniform across all outputs, since that is mathematically
%! ## equivalent to attributing the same delay to InputDelay instead.
%! ## Verify this succeeds (not "not yet supported") and matches the
%! ## InputDelay-based formulation exactly (ssdata/freqresp).
%! num = {1,1;1,1}; den = {[1 1],[1 2];[1 3],[1 4]};
%! sys_out = tf (num, den, "OutputDelay", [1.33;1.33]);
%! sys_in  = tf (num, den, "InputDelay", [1.33;1.33]);
%! d_out = c2d (sys_out, 0.5, "zoh");
%! d_in  = c2d (sys_in, 0.5, "zoh");
%! assert (isdt (d_out), true);
%! [Ao, Bo, Co, Do] = ssdata (d_out);
%! [Ai, Bi, Ci, Di] = ssdata (d_in);
%! assert (Ao, Ai, 1e-10);
%! assert (Bo, Bi, 1e-10);
%! assert (Co, Ci, 1e-10);
%! assert (Do, Di, 1e-10);
%! w = [0.1, 1, 5];
%! assert (freqresp (d_out, w), freqresp (d_in, w), 1e-10);

%!error <not yet supported> ...
%! ## row-varying OutputDelay (differs across outputs) is genuinely out of
%! ## scope -- it cannot be pushed through to the input side -- and must
%! ## still fall through to the "not yet supported" error.
%! c2d (tf ({1,1;1,1}, {[1 1],[1 2];[1 3],[1 4]}, "OutputDelay", [1.33;0.5]), 0.5, "zoh")

%!error <not yet supported> ...
%! ## row-varying IODelay within a column is likewise genuinely out of
%! ## scope and must still error.
%! sys = ss ([-1, 0; 0, -2], [1, 0; 0, 1], [1, 0; 0, 1], [0, 0; 0, 0]);
%! sys = set (sys, "IODelay", [1.33, 0; 0.5, 0]);
%! c2d (sys, 0.5, "zoh")


%!test  # c2dOptions struct as 3rd argument, default DelayModeling
%! opt = c2dOptions ();
%! sys = tf (1, [1 1], "InputDelay", 1.0);
%! dsys = c2d (sys, 0.5, opt);
%! assert (dsys.InputDelay, 2);

%!test  # DelayModeling='state': fractional SISO delay rounds to integer samples, absorbed into extra dynamics
%! opt = c2dOptions ("DelayModeling", "state");
%! sys = tf (1, [1 1], "InputDelay", 1.33);
%! dsys = c2d (sys, 0.5, opt);
%! assert (isdt (dsys), true);
%! assert (hasdelay (dsys), false);
%! rational_dsys = c2d (tf (1, [1 1]), 0.5, "zoh");
%! k = round (1.33 / 0.5);
%! w = [0.1, 1, 5];
%! expected = freqresp (rational_dsys, w) .* reshape (exp (-1i * w * 0.5 * k), 1, 1, []);
%! assert (freqresp (dsys, w), expected, 1e-8);

%!test  # DelayModeling='state': MIMO fractional delay succeeds (unlike 'delay' mode), absorbed per channel
%! opt = c2dOptions ("DelayModeling", "state");
%! sys = tf ({1,1;1,1}, {[1 1],[1 2];[1 3],[1 4]}, "InputDelay", [1.33;0]);
%! dsys = c2d (sys, 0.5, opt);
%! assert (hasdelay (dsys), false);
%! rational_dsys = c2d (tf ({1,1;1,1}, {[1 1],[1 2];[1 3],[1 4]}), 0.5, "zoh");
%! total = [round(1.33/0.5), 0; round(1.33/0.5), 0];
%! w = [0.1, 1, 5];
%! resp = freqresp (rational_dsys, w);
%! expected = zeros (2, 2, numel (w));
%! for i = 1:2
%!   for j = 1:2
%!     expected(i,j,:) = reshape (resp(i,j,:), 1, []) .* exp (-1i * w * 0.5 * total(i,j));
%!   endfor
%! endfor
%! assert (freqresp (dsys, w), expected, 1e-8);

%!test  # DelayModeling='state': ss input absorbs InputDelay into extra dynamics
%! sys = ss (-1, 1, 1, 0, "InputDelay", 0.65);
%! opt = c2dOptions ("DelayModeling", "state");
%! dsys = c2d (sys, 0.5, opt);
%! assert (isa (dsys, "ss"), true);
%! assert (isdt (dsys), true);
%! assert (hasdelay (dsys), false);
%! rational_dsys = c2d (ss (-1, 1, 1, 0), 0.5, "zoh");
%! k = round (0.65 / 0.5);
%! w = [0.1, 1, 5];
%! expected = freqresp (rational_dsys, w) .* reshape (exp (-1i * w * 0.5 * k), 1, 1, []);
%! assert (freqresp (dsys, w), expected, 1e-8);

%!test  # DelayModeling='state': ss input with per-channel IODelay matrix (MIMO)
%! A = [-1, 0; 0, -2];
%! B = [1, 0; 0, 1];
%! C = [1, 0; 0, 1];
%! D = [0, 0; 0, 0];
%! sys = ss (A, B, C, D);
%! sys = set (sys, "IODelay", [0.5, 0; 0, 1.0]);
%! opt = c2dOptions ("DelayModeling", "state");
%! dsys = c2d (sys, 0.5, opt);
%! assert (hasdelay (dsys), false);
%! rational_dsys = c2d (ss (A, B, C, D), 0.5, "zoh");
%! total = [round(0.5/0.5), 0; 0, round(1.0/0.5)];
%! w = [0.1, 1, 5];
%! resp = freqresp (rational_dsys, w);
%! expected = zeros (2, 2, numel (w));
%! for i = 1:2
%!   for j = 1:2
%!     expected(i,j,:) = reshape (resp(i,j,:), 1, []) .* exp (-1i * w * 0.5 * total(i,j));
%!   endfor
%! endfor
%! assert (freqresp (dsys, w), expected, 1e-8);

%!error <cannot combine> c2d (tf (1, [1 1]), 0.5, c2dOptions (), 100)

%!test  # exact-integer-sample InternalDelay: feedback()-produced system, tau/tsam is a whole number
%! T = 0.5;
%! G = ss (-1, 1, 1, 0, "IODelay", T);
%! L = feedback (G);
%! assert (hasinternaldelay (L), true);
%! dsys = c2d (L, 0.25, "zoh");
%! assert (isdt (dsys), true);
%! assert (hasinternaldelay (dsys), true);
%! assert (dsys.internaldelay, T / 0.25, 1e-10);

%!test  # fractional InternalDelay: rounds to nearest sample
%! T = 0.33;
%! tsam = 0.1;
%! G = ss (-1, 1, 1, 0, "IODelay", T);
%! L = feedback (G);
%! dsys = c2d (L, tsam, "zoh");
%! assert (hasinternaldelay (dsys), true);
%! assert (dsys.internaldelay, round (T / tsam));

%!test  # zoh sanity check: discretized InternalDelay system tracks continuous freqresp at low frequency
%! T = 0.3;
%! G = ss (-2, 1, 1, 0, "IODelay", T);
%! L = feedback (G);
%! tsam = 0.01;
%! dsys = c2d (L, tsam, "zoh");
%! w = [0.05, 0.2, 1];
%! Hc = freqresp (L, w);
%! Hd = freqresp (dsys, w);
%! assert (Hd, Hc, 5e-2);

%!test  # InternalDelay with no delay ports (b2/c2/d12/d21/d22 empty): degenerate extended system still discretizes
%! sys = set (ss (-1, 1, 1, 0), "internaldelay", 0.5);
%! dsys = c2d (sys, 0.5, "zoh");
%! assert (isdt (dsys), true);
%! assert (dsys.internaldelay, 1);

%!test  # combined InternalDelay AND ordinary InputDelay: both rounded to samples
%! T = 0.5;
%! tsam = 0.25;
%! G = ss (-1, 1, 1, 0, "IODelay", T);
%! L = feedback (G);
%! L = set (L, "InputDelay", 0.75);
%! assert (hasinternaldelay (L), true);
%! assert (hasdelay (L), true);
%! dsys = c2d (L, tsam, "zoh");
%! assert (isdt (dsys), true);
%! assert (hasinternaldelay (dsys), true);
%! assert (dsys.internaldelay, T / tsam, 1e-10);
%! assert (dsys.InputDelay, round (0.75 / tsam));


%!test  # MIMO InternalDelay: two channels, DIFFERENT delays, dimension-correct multi-port c2d
%! ## Genuine 2-in/2-out InternalDelay fixture built as append() of two
%! ## independent SISO feedback loops (the isolated-SISO-block topology
%! ## __sys_connect__ actually supports).  Distinct delays (0.4 vs 0.8) so an
%! ## index/port swap could not hide behind equal values.
%! T1 = 0.4; T2 = 0.8; tsam = 0.2;
%! G1 = ss (-1, 1, 1, 0, "IODelay", T1);
%! G2 = ss (-2, 1, 1, 0, "IODelay", T2);
%! sys = append (feedback (G1), feedback (G2));
%! assert (get (sys, "internaldelay"), [T1; T2], 1e-12);   # two channels
%! dsys = c2d (sys, tsam, "zoh");
%! assert (hasinternaldelay (dsys), true);
%! assert (isdt (dsys), true);
%! assert (get (dsys, "internaldelay"), [T1/tsam; T2/tsam], 1e-10);  # 2 and 4 samples
%! assert (size (dsys), [2, 2]);
%! ## Independent reference: each diagonal channel must equal the SISO c2d of
%! ## that loop alone; off-diagonals must be exactly zero (channels decoupled).
%! d1 = c2d (feedback (G1), tsam, "zoh");
%! d2 = c2d (feedback (G2), tsam, "zoh");
%! w = [0.1, 0.7, 2];
%! H = freqresp (dsys, w); H1 = freqresp (d1, w)(:); H2 = freqresp (d2, w)(:);
%! for k = 1:numel (w)
%!   assert (H(1,1,k), H1(k), 1e-12);
%!   assert (H(2,2,k), H2(k), 1e-12);
%!   assert (H(1,2,k), 0, 1e-12);
%!   assert (H(2,1,k), 0, 1e-12);
%! endfor

## InternalDelay that rounds to 0 samples must error clearly, not silently
## drop the port's feedthrough (see the tau>=1 assumption documented in
## __delay_lookup__/__buffered_sim__).
%!error <would round to 0 samples> c2d (feedback (ss (-1, 1, 1, 0, "IODelay", 0.01)), 1, "zoh")

%!test  # exact fractional IODelay via zoh: matches MATLAB's split-ZOH result
%! h = tf (10, [1 3 10], "IODelay", 0.25);
%! hd = c2d (h, 0.1);      # method defaults to "zoh"
%! [num, den] = tfdata (hd);
%! assert (hd.InputDelay + hd.IODelay, 3);   # z^-3 pure delay factored out
%! assert (num{1}, [0.01187, 0.06408, 0.009721], 1e-4);
%! assert (den{1}, [1, -1.655, 0.7408], 1e-3);

function tf_ok = __c2d_frac_applicable__ (sys, tsam)
  [outdelay] = get (sys, "outputdelay");
  total = totaldelay (sys);           # p-by-m matrix, seconds
  [p, m] = size (total);
  ## A single-output system's OutputDelay commutes with its InputDelay(s)
  ## (it is equivalent, for a scalar-output LTI system, to shifting every
  ## input by the same amount), so it can always be folded into the
  ## per-column totaldelay used below -- this is exactly what Task 1's
  ## SISO OutputDelay regression test relies on.  For p > 1 a nonzero,
  ## per-row-varying OutputDelay would make the per-column delay
  ## inconsistent across outputs; that case is already rejected by the
  ## uniformity check below, so no separate early-return is needed here.
  ##
  ## Deliberate scope decision (not just an accepted side effect): for
  ## MIMO (p > 1), OutputDelay may be nonzero as long as it is *uniform*
  ## across all outputs -- that is mathematically equivalent to
  ## attributing the same delay to InputDelay instead (every output being
  ## delayed by the same amount commutes with delaying every input by that
  ## amount), and is folded into the per-column totaldelay uniformity
  ## check below exactly like any other InputDelay/IODelay contribution.
  ## Only row-VARYING OutputDelay is rejected here, since that would
  ## require genuine output-side treatment (delaying different outputs by
  ## different amounts cannot be pushed through to the input side), which
  ## this fix does not implement -- such systems fall through to the
  ## "not yet supported" error below.
  if (p > 1 && any (outdelay != outdelay(1)))
    tf_ok = false;
    return;
  endif
  for j = 1 : m
    col = total (:, j);
    if (any (abs (col - col(1)) > 1e-10))
      tf_ok = false;
      return;
    endif
  endfor
  tf_ok = true;
endfunction

function sys = __c2d_frac_ss__ (origsys, tsam)
  sys_ss = ss (origsys);
  [A, B, C, D] = ssdata (sys_ss);
  total = totaldelay (origsys);
  [p, m] = size (total);
  n = rows (A);

  col_tau = total (1, :);       # one delay value per column (uniform down
                                 # each column, guaranteed by the
                                 # applicability check)

  extra_needed = false (1, m);
  d_samples = zeros (1, m);
  Bd_cols = cell (1, m);
  g0_cols = cell (1, m);

  for j = 1 : m
    [Ad, Bd, extra] = __c2d_frac_zoh__ (A, B(:,j), tsam, col_tau(j));
    d_samples(j) = extra.d;
    Bd_cols{j} = Bd;
    extra_needed(j) = extra.needed;
    if (extra.needed)
      g0_cols{j} = extra.g0;
      ## NOTE: no "+1" here.  The extra register state remains part of the
      ## augmented dynamics (A_aug/B_aug below) and its own eigenvalue-zero
      ## pole already contributes exactly one sample of delay to the
      ## transfer function; adding 1 to the explicit InputDelay field here
      ## as well would double-count that sample (verified against an
      ## independent Astrom & Wittenmark split-ZOH ground truth: with the
      ## register left in the dynamics, InputDelay must stay at d, not
      ## d+1).  The one case where the register's dynamics are explicitly
      ## removed from the transfer function (the SISO tf/zpk pole/zero
      ## cancellation below) re-adds the missing sample there, exactly
      ## because the register's contribution is being deleted from the
      ## dynamics at that point.
    endif
  endfor

  n_extra = sum (extra_needed);
  A_aug = zeros (n + n_extra, n + n_extra);
  A_aug (1:n, 1:n) = Ad;
  B_aug = zeros (n + n_extra, m);
  C_aug = [C, zeros(rows(C), n_extra)];

  row = n + 1;
  for j = 1 : m
    B_aug (1:n, j) = Bd_cols{j};
    if (extra_needed(j))
      A_aug (1:n, row) = g0_cols{j};
      B_aug (row, j) = 1;
      row += 1;
    endif
  endfor

  sys = ss (A_aug, B_aug, C_aug, D, tsam);
  sys = set (sys, "InputDelay", d_samples(:));

  if (isa (origsys, "tf") || isa (origsys, "zpk"))
    sys = tf (sys);

    ## The extra register state added above (when extra.needed) has an
    ## eigenvalue of exactly zero and contributes no real dynamics -- it
    ## exists purely to carry the sub-sample remainder through one more
    ## sample.  Converting the augmented state-space realization to a
    ## transfer function therefore leaves an exact pole/zero pair at the
    ## origin (a leading zero numerator coefficient paired with a trailing
    ## ~0 denominator coefficient): cancel it explicitly so num/den come
    ## back at minimal (MATLAB-compatible) order instead of one degree too
    ## high.
    if (p == 1 && m == 1 && n_extra > 0)
      [num, den] = tfdata (sys, "v");
      if (abs (num(1)) < 1e-8 * max (abs (num)) ...
          && abs (den(end)) < 1e-8 * max (abs (den)))
        num = num(2:end);
        den = den(1:end-1);
        sys = tf (num, den, tsam);
        ## The register's dynamics were just deleted from num/den above, so
        ## (unlike the general case) the one sample of delay it used to
        ## contribute via its own pole must now be re-added explicitly.
        sys = set (sys, "InputDelay", d_samples + 1);
      endif
    endif

    if (isa (origsys, "zpk"))
      sys = zpk (sys);
    endif
  endif
endfunction
