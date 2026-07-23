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
## Frequency response of SS models.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.6

function H = __freqresp__ (sys, w, cellflag = false)

  ## Do not prescale internal-delay systems: prescale() rebuilds the ss without
  ## the delay ports (b2/c2/d12/d21/d22/internaldelay), which would silently drop
  ## the internal delay.  Prescaling is only a numerical-conditioning nicety.
  if (sys.scaled == false && ! hasinternaldelay (sys))
    sys = prescale (sys);
  endif

  [a, b, c, d, e, tsam] = dssdata (sys);

  if (isct (sys))  # continuous system
    s = i * w;
  else             # discrete system
    s = exp (i * w * abs (tsam));
  endif

  ## compute H and prevent warnings when det(e*jw-A) gets close
  ## to a singular matrix near to resonance frequencies
  warning_id_near_sing = "Octave:nearly-singular-matrix";
  warning_near_sing = warning ("query", warning_id_near_sing);
  warning ("off", warning_id_near_sing);

  if (hasinternaldelay (sys))
    ## LFT closure for internal delay: the ordinary dynamics (a,b,c,d) are the
    ## (A,B1,C1,D11) block; the delay ports are (b2,c2,d12,d21,d22).  The delay
    ## law is w = Delta z with Delta = diag(exp(-jw*tau)) (continuous) or
    ## diag(exp(-jw*tsam*tau)) (discrete), tau = sys.internaldelay.
    ##   T11 = D11 + C1 (jwI-A)^-1 B1     T12 = D12 + C1 (jwI-A)^-1 B2
    ##   T21 = D21 + C2 (jwI-A)^-1 B1     T22 = D22 + C2 (jwI-A)^-1 B2
    ##   H   = T11 + T12 Delta (I - T22 Delta)^-1 T21
    b2  = sys.b2;
    c2  = sys.c2;
    d12 = sys.d12;
    d21 = sys.d21;
    d22 = sys.d22;
    tau = sys.internaldelay(:);
    nd  = numel (tau);
    Ind = eye (nd);
    if (isct (sys))
      delta = @(wk) exp (-1i * wk * tau);
    else
      delta = @(wk) exp (-1i * wk * abs (tsam) * tau);
    endif
    H = cell (size (s));
    for k = 1 : numel (s)
      M   = s(k)*e - a;
      T11 = c /M*b  + d;
      T12 = c /M*b2 + d12;
      T21 = c2/M*b  + d21;
      T22 = c2/M*b2 + d22;
      D   = diag (delta (w(k)));
      H{k} = T11 + T12*D*((Ind - T22*D) \ T21);
    endfor
  else
    H = arrayfun (@(x) c/(x*e - a)*b + d, s, "uniformoutput", false);
  endif

  warning (warning_near_sing.state, warning_id_near_sing);

  if (! cellflag)
    H = cat (3, H{:});
  endif

endfunction


%!test  # no internal delay: byte-for-byte the ordinary D+C(jwI-A)^-1 B path
%! sys = ss (-2, 1, 1, 0);
%! w = [0.1, 1, 5, 20];
%! H = __freqresp__ (sys, w);
%! a = -2; b = 1; c = 1; d = 0;
%! ref = arrayfun (@(x) c/(1i*x - a)*b + d, w, "uniformoutput", false);
%! ref = cat (3, ref{:});
%! assert (H, ref);

%!test  # scalar unity-feedback loop with IODelay: 1/(1+G e^{-jwT})-style closure
%! ## G(s) = 1/(s+1) with a loop delay T = 0.3 s; L = feedback (G) closes
%! ## the loop, trapping the delay internally.  Closed loop:
%! ##   H(jw) = G e^{-jwT} / (1 + G e^{-jwT}),  G = 1/(jw+1).
%! T = 0.3;
%! G = ss (-1, 1, 1, 0, "IODelay", T);
%! L = feedback (G);
%! assert (hasinternaldelay (L));
%! assert (get (L, "internaldelay"), T, 1e-12);
%! w = [0.05, 0.5, 1.7, 4, 11];
%! H = __freqresp__ (L, w);
%! H = reshape (H, 1, []);
%! Gd = (1 ./ (1i*w + 1)) .* exp (-1i * w * T);
%! expected = Gd ./ (1 + Gd);
%! assert (H, expected, 1e-9);

%!test  # cellflag=true returns the same values in cell form
%! T = 0.25;
%! G = ss (-1, 1, 1, 0, "IODelay", T);
%! L = feedback (G);
%! w = [0.3, 2, 6];
%! Hc = __freqresp__ (L, w, true);
%! Ha = __freqresp__ (L, w, false);
%! for k = 1:numel (w)
%!   assert (Hc{k}, Ha(:,:,k), 1e-12);
%! endfor

%!test  # freqresp() top-level wrapper dispatches to __freqresp__ on an internal-delay system
%! T = 0.4;
%! G = ss (-2, 1, 1, 0, "IODelay", T);
%! L = feedback (G);
%! w = [0.1, 1, 3];
%! assert (freqresp (L, w), __freqresp__ (L, w), 1e-12);

%!test  # MIMO InternalDelay: two channels with DIFFERENT delays, LFT closure per channel
%! ## Genuine 2-in/2-out fixture (append of two independent SISO feedback loops)
%! ## so the internal-delay LFT closure is exercised for more than one port at
%! ## once.  Reference is the hand-derived closed form of each decoupled channel
%! ## H_k(jw) = Gk e^{-jwTk} / (1 + Gk e^{-jwTk}),  Gk = 1/(jw+ak); off-diagonals
%! ## must be identically zero.
%! T1 = 0.3; a1 = 1; T2 = 0.7; a2 = 2;
%! G1 = ss (-a1, 1, 1, 0, "IODelay", T1);
%! G2 = ss (-a2, 1, 1, 0, "IODelay", T2);
%! sys = append (feedback (G1), feedback (G2));
%! assert (get (sys, "internaldelay"), [T1; T2], 1e-12);
%! w = [0.05, 0.5, 1.7, 4, 11];
%! H = __freqresp__ (sys, w);
%! Gd1 = (1 ./ (1i*w + a1)) .* exp (-1i*w*T1); e1 = Gd1 ./ (1 + Gd1);
%! Gd2 = (1 ./ (1i*w + a2)) .* exp (-1i*w*T2); e2 = Gd2 ./ (1 + Gd2);
%! for k = 1:numel (w)
%!   assert (H(1,1,k), e1(k), 1e-10);
%!   assert (H(2,2,k), e2(k), 1e-10);
%!   assert (H(1,2,k), 0, 1e-12);
%!   assert (H(2,1,k), 0, 1e-12);
%! endfor
