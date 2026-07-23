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
## @deftypefn {Function File} {@var{sys} =} feedback (@var{sys1})
## @deftypefnx {Function File} {@var{sys} =} feedback (@var{sys1}, @var{"+"})
## @deftypefnx {Function File} {@var{sys} =} feedback (@var{sys1}, @var{sys2})
## @deftypefnx {Function File} {@var{sys} =} feedback (@var{sys1}, @var{sys2}, @var{"+"})
## @deftypefnx {Function File} {@var{sys} =} feedback (@var{sys1}, @var{sys2}, @var{feedin}, @var{feedout})
## @deftypefnx {Function File} {@var{sys} =} feedback (@var{sys1}, @var{sys2}, @var{feedin}, @var{feedout}, @var{"+"})
## Feedback connection of two @acronym{LTI} models.
##
## @strong{Inputs}
## @table @var
## @item sys1
## @acronym{LTI} model of forward transmission.  @code{[p1, m1] = size (sys1)}.
## @item sys2
## @acronym{LTI} model of backward transmission.
## If not specified, an identity matrix of appropriate size is taken.
## @item feedin
## Vector containing indices of inputs to @var{sys1} which are involved in the feedback loop.
## The number of @var{feedin} indices and outputs of @var{sys2} must be equal.
## If not specified, @code{1:m1} is taken.
## @item feedout
## Vector containing indices of outputs from @var{sys1} which are to be connected to @var{sys2}.
## The number of @var{feedout} indices and inputs of @var{sys2} must be equal.
## If not specified, @code{1:p1} is taken.
## @item "+"
## Positive feedback sign.  If not specified, @var{"-"} for a negative feedback interconnection
## is assumed.  @var{+1} and @var{-1} are possible as well, but only from the third argument
## onward due to ambiguity.
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## Resulting @acronym{LTI} model.
## @end table
##
## @strong{Remarks}
## Both operands may contain input, output, and I/O delays.  If a delay is trapped
## inside the closed loop, it will appear as an @code{InternalDelay} property on the
## (always-@code{ss}) result.
##
## @strong{Block Diagram}
## @example
## @group
##  u    +         +--------+             y
## ------>(+)----->|  sys1  |-------+------->
##         ^ -     +--------+       |
##         |                        |
##         |       +--------+       |
##         +-------|  sys2  |<------+
##                 +--------+
## @end group
## @end example
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.7

function sys = feedback (sys1, sys2, feedin, feedout, fbsign = -1)

  [p1, m1] = size (sys1);

  switch (nargin)
    case 1                          # sys = feedback (sys)
      if (p1 != m1)
        error ("feedback: argument must be a square system");
      endif
      sys2 = eye (p1);
      feedin = 1 : m1;
      feedout = 1 : p1;

    case 2
      if (ischar (sys2))            # sys = feedback (sys, "+")
        if (p1 != m1)
          error ("feedback: first argument must be a square system");
        endif
        fbsign = __check_fbsign__ (sys2);
        sys2 = eye (p1);
      endif                         # sys = feedback (sys1, sys2)
      feedin = 1 : m1;
      feedout = 1 : p1;

    case 3                          # sys = feedback (sys1, sys2, "+")
      fbsign = __check_fbsign__ (feedin);
      feedin = 1 : m1;
      feedout = 1 : p1;

    case 4                          # sys = feedback (sys1, sys2, feedin, feedout)
      ## nothing needs to be done here
      ## case 4 required to prevent "otherwise"

    case 5                          # sys = feedback (sys1, sys2, feedin, feedout, "+")
      fbsign = __check_fbsign__ (fbsign);

    otherwise
      print_usage ();
  endswitch

  if (ischar (feedin))
    feedin = {feedin};
  endif
  if (ischar (feedout))
    feedout = {feedout};
  endif

  if (iscell (feedin))
    tmp = cellfun (@(x) __str2idx__ (sys1.ingroup, sys1.inname, x, "in"), feedin, "uniformoutput", false);
    feedin = vertcat (tmp{:});
  endif
  if (iscell (feedout))
    tmp = cellfun (@(x) __str2idx__ (sys1.outgroup, sys1.outname, x, "out"), feedout, "uniformoutput", false);
    feedout = vertcat (tmp{:});
  endif

  if (! is_real_vector (feedin) || ! isequal (feedin, abs (fix (feedin))))
    error ("feedback: require 'feedin' to be a vector of integers");
  endif
  if (! is_real_vector (feedout) || ! isequal (feedout, abs (fix (feedout))))
    error ("feedback: require 'feedout' to be a vector of integers");
  endif

  [p2, m2] = size (sys2);

  l_feedin = length (feedin);
  l_feedout = length (feedout);

  if (l_feedin != p2)
    error ("feedback: feedin indices: %d, outputs sys2: %d", l_feedin, p2);
  endif
  if (l_feedout != m2)
    error ("feedback: feedout indices: %d, inputs sys2: %d", l_feedout, m2);
  endif

  if (any (feedin > m1 | feedin < 1))
    error ("feedback: range of feedin indices exceeds dimensions of sys1");
  endif
  if (any (feedout > p1 | feedout < 1))
    error ("feedback: range of feedout indices exceeds dimensions of sys1");
  endif

  M11 = zeros (m1, p1);
  M22 = zeros (m2, p2);

  M12 = full (sparse (feedin, 1:l_feedin, fbsign, m1, p2));
  M21 = full (sparse (1:l_feedout, feedout, 1, m2, p1));
  
  ## NOTE: for-loops do NOT the same as
  ##       M12(feedin, 1:l_feedin) = fbsign;
  ##       M21(1:l_feedout, feedout) = 1;
  ##
  ## M12 = zeros (m1, p2);
  ## M21 = zeros (m2, p1);
  ##
  ## for k = 1 : l_feedin
  ##   M12(feedin(k), k) = fbsign;
  ## endfor
  ##
  ## for k = 1 : l_feedout
  ##   M21(k, feedout(k)) = 1;
  ## endfor

  M = [M11, M12;
       M21, M22];

  in_idx = 1 : m1;
  out_idx = 1 : p1;

  ## A delay trapped inside the loop cannot be represented by tf/zpk, which
  ## carry no InternalDelay.  Force both operands to state-space so the delay
  ## survives as an internal delay port (matches Matlab: such a result is ss).
  d1 = isa (sys1, "lti") && (hasdelay (sys1) || hasinternaldelay (sys1));
  d2 = isa (sys2, "lti") && (hasdelay (sys2) || hasinternaldelay (sys2));
  if (d1 || d2)
    sys1 = ss (sys1);
    sys2 = ss (sys2);
  endif

  sys = __sys_group__ (sys1, sys2);
  sys = __sys_connect__ (sys, M);
  sys = __sys_prune__ (sys, out_idx, in_idx);

endfunction


function fbsign = __check_fbsign__ (fbsign)

  if (is_real_scalar (fbsign))
    fbsign = sign (fbsign);
  elseif (ischar (fbsign))
    if (strcmp (fbsign, "+"))
      fbsign = +1;
    elseif (strcmp (fbsign, "-"))
      fbsign = -1;
    else
      error ("feedback: invalid feedback sign string");
    endif
  else
    error ("feedback: invalid feedback sign type");
  endif

endfunction


## Feedback inter-connection of two systems in state-space form
## Test from SLICOT AB05ND
%!shared M, Me
%! A1 = [ 1.0   0.0  -1.0
%!        0.0  -1.0   1.0
%!        1.0   1.0   2.0 ];
%!
%! B1 = [ 1.0   1.0   0.0
%!        2.0   0.0   1.0 ].';
%!
%! C1 = [ 3.0  -2.0   1.0
%!        0.0   1.0   0.0 ];
%!
%! D1 = [ 1.0   0.0
%!        0.0   1.0 ];
%!
%! A2 = [-3.0   0.0   0.0
%!        1.0   0.0   1.0
%!        0.0  -1.0   2.0 ];
%!
%! B2 = [ 0.0  -1.0   0.0
%!        1.0   0.0   2.0 ].';
%!
%! C2 = [ 1.0   1.0   0.0
%!        1.0   1.0  -1.0 ];
%!
%! D2 = [ 1.0   1.0
%!        0.0   1.0 ];
%!
%! sys1 = ss (A1, B1, C1, D1);
%! sys2 = ss (A2, B2, C2, D2);
%! sys = feedback (sys1, sys2);
%! [A, B, C, D] = ssdata (sys);
%! M = [A, B; C, D];
%!
%! Ae = [-0.5000  -0.2500  -1.5000  -1.2500  -1.2500   0.7500
%!       -1.5000  -0.2500   0.5000  -0.2500  -0.2500  -0.2500
%!        1.0000   0.5000   2.0000  -0.5000  -0.5000   0.5000
%!        0.0000   0.5000   0.0000  -3.5000  -0.5000   0.5000
%!       -1.5000   1.2500  -0.5000   1.2500   0.2500   1.2500
%!        0.0000   1.0000   0.0000  -1.0000  -2.0000   3.0000 ];
%!
%! Be = [ 0.5000   0.7500
%!        0.5000  -0.2500
%!        0.0000   0.5000
%!        0.0000   0.5000
%!       -0.5000   0.2500
%!        0.0000   1.0000 ];
%!
%! Ce = [ 1.5000  -1.2500   0.5000  -0.2500  -0.2500  -0.2500
%!        0.0000   0.5000   0.0000  -0.5000  -0.5000   0.5000 ];
%!
%! De = [ 0.5000  -0.2500
%!        0.0000   0.5000 ];
%!
%! Me = [Ae, Be; Ce, De];
%!
%!assert (M, Me, 1e-4);


## sensitivity function
%!shared S1, S2
%! P = ss (-2, 3, 4, 5);  # meaningless numbers
%! C = ss (-1, 1, 1, 0);  # ditto
%! L = P * C;
%! I = eye (size (L));
%! S1 = feedback (I, L);
%! S2 = inv (I + L);
%!assert (S1.a, S2.a, 1e-4);
%!assert (S1.b, S2.b, 1e-4);
%!assert (S1.c, S2.c, 1e-4);
%!assert (S1.d, S2.d, 1e-4);


%!test  # no delay on either operand: existing behavior unaffected (regression)
%! s1 = tf (1, [1 1]);
%! s = feedback (s1);
%! assert (hasdelay (s), false);

## Internal delay: an input delay closing a feedback loop is no longer an
## error; it is absorbed into the closed loop's InternalDelay (the delay
## e^{-0.1 s} lives inside the loop 1/(s+1+e^{-0.1 s}) and cannot be
## expressed as a pure I/O time-shift).
%!test
%! s = feedback (tf (1, [1 1], "InputDelay", 0.1));
%! assert (isa (s, "ss"));                       # forced to ss to carry delay
%! assert (get (s, "internaldelay"), 0.1, 1e-10);
%! assert (hasinternaldelay (s), true);

## Scalar loop with an IODelay in the forward path.  A single scalar delay
## closing a loop is irreducible, so the closed-loop InternalDelay is just T.
## Hand derivation:  G = e^{-sT} g(s),  C = c(s),  L = G C.
## Closed loop feedback(G,C) = G/(1+GC); the e^{-sT} sits inside the loop,
## so InternalDelay = T (0.3 here) and no residual I/O delay remains.
%!test
%! T = 0.3;
%! G = ss (-1, 1, 1, 0, "IODelay", T);
%! C = ss (-2, 1, 1, 0);
%! cl = feedback (G, C);
%! assert (get (cl, "internaldelay"), T, 1e-10);
%! assert (get (cl, "iodelay"), 0, 1e-10);       # absorbed, no residual

## feedback() and connect() must agree on InternalDelay for equivalent
## topologies, since they share __sys_connect__.
%!test
%! T = 0.3;
%! G = ss (-1, 1, 1, 0, "IODelay", T);
%! C = ss (-2, 1, 1, 0);
%! cl_fb = feedback (G, C);
%! ## build the same negative feedback loop with connect():
%! ##   inputs  : [r ; u_G ; u_C] = grouped [G.in ; C.in]  plus external r
%! ##   here use index-based cm on append (G, C):
%! ##   in 1 -> G, in 2 -> C ; out 1 -> G, out 2 -> C
%! Gall = append (G, C);
%! ## u_G = r - y_C  (row: input 1 fed by -output 2)
%! ## u_C = y_G      (row: input 2 fed by +output 1)
%! cm = [1, -2; 2, 1];
%! cl_cn = connect (Gall, cm, 1, 1);
%! assert (get (cl_cn, "internaldelay"), get (cl_fb, "internaldelay"), 1e-10);

## MIMO plant whose IODelay column has more than one differently-valued
## nonzero entry (one input feeding two differently-delayed outputs).  A single
## shared actuator column cannot carry two different delays, so the column is
## decomposed onto an independent shadow-state block (one per decomposed
## column) and each output taps it with its own delay.  This used to error
## ("not yet supported"); it is now realised exactly.  Coupled A so both
## outputs genuinely depend on the shared input.
%!test
%! a = [-1, 0.4; 0, -2];
%! b = [1; 1];                     # single input drives both states
%! c = eye (2);
%! d = [0; 0];
%! iod = [0.2; 0.5];               # one input, two differently-delayed outputs
%! G = ss (a, b, c, d, "IODelay", iod);
%! cl = feedback (G, [1, 1]);      # close the single loop (sum of outputs)
%! assert (rows (ssdata (cl)), 4); # 2 original + 2 shadow (one decomposed col)
%! assert (sort (get (cl, "internaldelay")), [0.2; 0.5], 1e-12);
%! w = [0.5, 1.3, 3.0];
%! for k = 1 : numel (w)
%!   s = 1i*w(k);
%!   T = c/(s*eye(2) - a)*b + d;            # 2x1 open-loop transfer
%!   Gt = T .* exp (-1i*w(k)*iod);          # per-entry delayed
%!   Href = Gt / (1 + [1, 1]*Gt);           # feedback (G, [1 1]): e = u - [1 1]y
%!   assert (freqresp (cl, w(k)), Href, 1e-9);
%! endfor
