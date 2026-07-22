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

  if (hasinternaldelay (sys))
    error ("c2d: InternalDelay is not yet supported");
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

      if (p != 1 || m != 1)
        error ("c2d: fractional delays on MIMO systems are not yet supported (per-channel Thiran approximation is not yet implemented)");
      endif

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


%!test  # SISO tf fractional InputDelay: approximated via thiran
%! sys = tf (1, [1 1], "InputDelay", 1.33);
%! dsys = c2d (sys, 0.5, "zoh");
%! assert (isdt (dsys), true);
%! assert (hasdelay (dsys), false);
%! sys_dyn = c2d (tf (1, [1 1]), 0.5, "zoh");
%! filt = thiran (1.33, 0.5);
%! expected = sys_dyn * filt;
%! [numd, dend] = tfdata (dsys);
%! [nume, dene] = tfdata (expected);
%! assert (numd, nume, 1e-8);
%! assert (dend, dene, 1e-8);

%!test  # SISO zpk fractional OutputDelay: approximated via thiran
%! sys = zpk ([], -1, 1, "OutputDelay", 1.33);
%! dsys = c2d (sys, 0.5, "zoh");
%! assert (isdt (dsys), true);
%! assert (hasdelay (dsys), false);
%! sys_dyn = c2d (zpk ([], -1, 1), 0.5, "zoh");
%! filt = thiran (1.33, 0.5);
%! expected = sys_dyn * filt;
%! w = [0.1, 1, 5];
%! assert (freqresp (dsys, w), freqresp (expected, w), 1e-8);

%!error <MIMO> c2d (tf ({1,1;1,1}, {[1 1],[1 2];[1 3],[1 4]}, "InputDelay", [1.33;0]), 0.5, "zoh")


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

%!error <InternalDelay> c2d (set (ss (-1, 1, 1, 0), "internaldelay", 0.5), 0.5, "zoh")

