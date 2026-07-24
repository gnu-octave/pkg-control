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

## Author: Prateek Ganguli <prateek.ganguli@gmail.com>
## Created: July 2026
## Version: 0.1

function H = __apply_freqresp_delay__ (sys, H, w, cellflag)

  if (! hasdelay (sys))
    return;
  endif

  total = totaldelay (sys);

  if (isdt (sys))
    tsam = abs (get (sys, "tsam"));
    phase = @(wk) exp (-1i * wk * tsam * total);
  else
    phase = @(wk) exp (-1i * wk * total);
  endif

  if (cellflag)
    for k = 1 : numel (w)
      H{k} = H{k} .* phase (w(k));
    endfor
  else
    for k = 1 : numel (w)
      H(:,:,k) = H(:,:,k) .* phase (w(k));
    endfor
  endif

endfunction


%!test  # no-delay system: exact no-op, cellflag=false
%! sys = tf (1, [1, 1]);
%! w = [0.1, 1, 5];
%! H = __freqresp__ (sys, w);
%! H2 = __apply_freqresp_delay__ (sys, H, w, false);
%! assert (H2, H);

%!test  # no-delay system: exact no-op, cellflag=true
%! sys = tf (1, [1, 1]);
%! w = [0.1, 1, 5];
%! H = __freqresp__ (sys, w, true);
%! H2 = __apply_freqresp_delay__ (sys, H, w, true);
%! assert (H2, H);

%!test  # SISO continuous tf with InputDelay, cellflag=false
%! sys = tf (1, [1, 1], "InputDelay", 0.3);
%! rational = tf (1, [1, 1]);
%! w = [0.1, 1, 5];
%! H = __freqresp__ (sys, w);
%! H2 = __apply_freqresp_delay__ (sys, H, w, false);
%! expected = __freqresp__ (rational, w) .* reshape (exp (-1i * w * 0.3), 1, 1, []);
%! assert (H2, expected, 1e-10);

%!test  # SISO discrete tf with integer InputDelay (samples), cellflag=false
%! sys = tf (1, [1, -0.5], 0.1, "InputDelay", 2);
%! rational = tf (1, [1, -0.5], 0.1);
%! w = [0.1, 1, 5];
%! H = __freqresp__ (sys, w);
%! H2 = __apply_freqresp_delay__ (sys, H, w, false);
%! expected = __freqresp__ (rational, w) .* reshape (exp (-1i * w * 0.1 * 2), 1, 1, []);
%! assert (H2, expected, 1e-10);

%!test  # MIMO tf with per-channel IODelay matrix, cellflag=false
%! sys = tf ({1, 1; 1, 1}, {[1, -0.5], [1, -0.6]; [1, -0.7], [1, -0.8]}, 0.1);
%! sys = set (sys, "IODelay", [1, 0; 0, 2]);
%! rational = tf ({1, 1; 1, 1}, {[1, -0.5], [1, -0.6]; [1, -0.7], [1, -0.8]}, 0.1);
%! total = [1, 0; 0, 2];
%! w = [0.1, 1, 5];
%! H = __freqresp__ (sys, w);
%! H2 = __apply_freqresp_delay__ (sys, H, w, false);
%! resp = __freqresp__ (rational, w);
%! expected = zeros (2, 2, numel (w));
%! for i = 1:2
%!   for j = 1:2
%!     expected(i,j,:) = reshape (resp(i,j,:), 1, []) .* exp (-1i * w * 0.1 * total(i,j));
%!   endfor
%! endfor
%! assert (H2, expected, 1e-8);

%!test  # cellflag=true with delay: SISO continuous tf
%! sys = tf (1, [1, 1], "InputDelay", 0.3);
%! rational = tf (1, [1, 1]);
%! w = [0.1, 1, 5];
%! H = __freqresp__ (sys, w, true);
%! H2 = __apply_freqresp_delay__ (sys, H, w, true);
%! rational_H = __freqresp__ (rational, w, true);
%! for k = 1:numel (w)
%!   expected_k = rational_H{k} * exp (-1i * w(k) * 0.3);
%!   assert (H2{k}, expected_k, 1e-10);
%! endfor
