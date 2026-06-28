## Copyright (C) 2025-2026 Torten Lilge
##
## This function is part of the GNU Octave Control Package
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {@var{info} =} stepinfo (@var{sys})
##
## Calculate some characteristic values of a system's step response
##
## @strong{Inputs}
## @table @var
## @item sys
## LTI system
## @end table
##
## @strong{Outputs}
## @table @var
## @item info
## Structure with the following fields
## @end table
##
## @seealso{step}
## @end deftypefn

function info = stepinfo (varargin)

  ## Check input arguments

  if (nargin < 1)
    error ("stepinfo: at least one input argument required\n");
  endif

  sys = varargin{1};
  if (! isa (sys, "lti"))
    error ("stepinfo: first argument must be an lti system\n");
  endif

  ## Thresholds

  thresholds = struct ();
  thresholds.t_rise = [0.1 0.9];
  thresholds.t_settling = 0.02;
  thresholds.t_transient = 0.02;

  ## Get information on the system

  [p, m] = size (sys);
  [num, den] = tfdata (sys);
  tfcell = cellfun (@(n,d) tf(n,d,sys.ts), num, den, "uniformoutput", false);

  stability = mat2cell (cellfun (@isstable, tfcell), ...
                        ones (1,p),  ones (1,m));
  ct = isct (sys);

  [y,t] = step (sys);
  N = length (t);

  y_cell = mat2cell (reshape (y,p*N,m), N*ones(1,p), ones (1,m));

  K = dcgain (sys);
  y_final = K;
  y_init = 0 * y_final; # for now, only y_init = 0 is considered
  y_final = mat2cell (y_final, ones(1,p), ones (1,m));
  y_init = mat2cell (y_init, ones(1,p), ones (1,m));

  ## RiseTime
  t_rise = cellfun (@(yc,yf,yi,st) __get_rise_time__ (yc,yf,yi,st,ct,t,thresholds.t_rise), ...
                    y_cell, y_final, y_init, stability, "uniformoutput", false);
  info = struct ("RiseTime", ...
                 cellfun (@diff, t_rise, "uniformoutput", false));

  ## TransientTime
  t_transient = cellfun (@(yc,yf,yi,st) __get_transient_time__ (yc,yf,yi,st,ct,t,'t',thresholds.t_transient), ...
                         y_cell, y_final, y_init, stability, "uniformoutput", false);
  info = __add_field__ (info, "TransientTime", t_transient);

  ## SettlingTime
  t_settling = cellfun (@(yc,yf,yi,st) __get_transient_time__ (yc,yf,yi,st,ct,t,'s',thresholds.t_settling), ...
                        y_cell, y_final, y_init, stability, "uniformoutput", false);
  info = __add_field__ (info, "SettlingTime", t_settling);

  ## SettlingMin/Max
  risen_idx = cellfun (@(tr) find(t >= tr(2))(1), t_rise, "uniformoutput", false);
  t_risen = cellfun (@(idx) t(idx:end), risen_idx, "uniformoutput", false);
  y_risen = cellfun (@(yc,idx) yc(idx:end), y_cell, risen_idx, ...
                     "uniformoutput", false);
  if (ct)
    y_risen = cellfun (@(yr) [thresholds.t_rise(2)*y_final; yr], y_risen, "uniformoutput", false);
  endif

  info = __add_field__ (info, "SettlingMin", cellfun (@min, y_risen, ...
                                             "uniformoutput", false));
  info = __add_field__ (info, "SettlingMax", cellfun (@max, y_risen, ...
                                             "uniformoutput", false));

  ## Overshoot / Undershoot
  y_rel = cellfun (@(yc,yi) (yc-yi), y_cell, y_init, "uniformoutput", false);
  y_norm = cellfun (@(yr,yi,yf) yr./(yf-yi), y_rel, y_init, y_final, ...
                    "uniformoutput", false);

  info = __add_field__ (info, "Overshoot", cellfun (@(yn) max(0,100*max(yn-1)), y_norm, ...
                                           "uniformoutput", false));
  info = __add_field__ (info, "Undershoot", cellfun (@(yn) min(0,-100*min(yn)), y_norm, ...
                                            "uniformoutput", false));

  ## Peak and Peak Time
  y_peak = cellfun (@(yr,yi,yf) sign(yf-yi)*yr, y_rel, y_init, y_final, ...
                    "uniformoutput", false);
  peak = cellfun (@max, y_peak, "uniformoutput", false);
  info = __add_field__ (info, "Peak", peak);

  peak_time = cellfun (@(yp,p) t(find(yp==p)(1)), y_peak, peak, ...
                       "uniformoutput", false);
  info = __add_field__ (info, "PeakTime", peak_time);

endfunction



function t_rise = __get_rise_time__ (y, y_final, y_init, stability, ct, t, t_tol)

  if (! stability)
    return;
  endif

  y_rise = arrayfun (@(x) t_tol * x, y_final, "uniformoutput", false);
  t_rise = [-1 -2];

  s = sign (y_final - y_init);

  for i = 1:2
    j = (find (s*y > s*y_rise{1}(i)))(1);

    if (! ct)
      t_rise(i) = t(j);
      continue
    endif

    if (j > 1)
      t_rise(i) = interp1 ([y(j-1) y(j)], [t(j-1) t(j)], y_rise{1}(i));
    else
      if (i == 1)
        t_rise(i) = 0;
      else
        t_rise(i) = Inf;
      endif
    endif

  endfor

endfunction



function t_transient = __get_transient_time__ (y, y_final, y_init, stability, ct, t, t_type, t_tol)

  if (! stability)
    return;
  endif

  y_error = abs (y - y_final);
  if (t_type == 't')
    ## TransientTime
    y_tol = t_tol * max (y_error);
  else
    ## SettlingTime
    y_tol = t_tol * abs (y_final - y_init);
  endif

  idx = find (y_error > y_tol);

  if length (idx) > 0

    idx = idx(end);

    if (idx < length (y_error))

      if (! ct)
        t_transient = t(idx+1);
        return;
      endif

      t_transient = interp1 ([y_error(idx) y_error(idx+1)], [t(idx) t(idx+1)], y_tol);

    else

      t_transient = Inf;

    endif

  else

    t_transient = Inf;

  endif

endfunction




function s_new = __add_field__ (s_old, field, value)

  s_new = arrayfun(@(s,k) setfield(s, field, value{k}), ...
                   s_old(:), (1:numel(s_old))', 'UniformOutput', true);
  s_new = reshape(s_new, size (value));

endfunction

