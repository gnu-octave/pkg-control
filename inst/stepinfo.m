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
## @deftypefnx {Function File} {@var{info} =} stepinfo (@var{sys, @var{property}, @var{value}, ...  })
##
## Calculate some characteristic values of a system's step response
##
## The following symbols are used for explaning the characteristic
## values in the following help text:
##
## @itemize @bullet{}
## @item
## @tex
## \(y_i\)
## @end tex
## @ifnottex
## @code{yi}
## @end ifnottex
## : Initial output value at
## @tex
## \(t=0^-\),
## @end tex
## @ifnottex
## @code{t=-0},
## @end ifnottex
## before the step is applied.
## @item
## @tex
## \(y_f\)
## @end tex
## @ifnottex
## @code{yf}
## @end ifnottex
## : Final steady state output value.
## @item
## @tex
## \(y_e\)
## @end tex
## @ifnottex
## @code{ye}
## @end ifnottex
## : Absolute error between output
## @tex
## \(y\)
## @end tex
## @ifnottex
## @code{y}
## @end ifnottex
##  and final value
## @tex
## \(y_f\).
## @end tex
## @ifnottex
## @code{yf}.
## @end ifnottex
## @end itemize
##
## @strong{Inputs}
##
## @table @var
## @item sys
## LTI system
## @item property
## @item value
## Properties/value pairs changing some default thresholds.
## @table @samp
## @item RaiseTimeLimits
## Array with two entries with relaitve values of @inlinefmtifelse{tex, @code{y_f}, @code{yf}}, between which the
## rise time is determined. The values have to be between 0 and 1. The default is
## @tex
## \([r_L \, r_U] = [0.1 \, 0.9]\)
## @end tex
## @ifnottex
## @code{[rL rU] = [0.1 0.9]}
## @end ifnottex
## if this property is omitted.
## @item SettlingTimeThreshold
## Relative value of the absolute maximum or initial difference to @inlinefmtifelse{tex, @code{y_f}, @code{yf}}
## defining a tolerance range around @inlinefmtifelse{tex, @code{y_f}, @code{yf}} for TransientTime or
## SettlingTime (see below). The default is
## @tex
## \(s_T = 0.02\)
## @end tex
## @ifnottex
## @code{sT = 0.02}
## @end ifnottex
## if this property is omitted.
## @end table
## @end table
##
## @strong{Outputs}
##
## @table @var
## @item info
## Structure (or structure array for MIMO systems)
## with some characteristics of the step response of
## system @var{sys}. In particular, the fields in @var{info} are:
## @table @samp
## @item RaiseTime
## Time required from
## @tex
## \(r_L y_f\)
## @end tex
## @ifnottex
## @code{rL*yf}
## @end ifnottex
## to
## @tex
## \(r_U y_f\).
## @end tex
## @ifnottex
## @code{rU*yf}.
## @end ifnottex
## @item TransientTime
## Time after which the absolute error
## @tex
## \(y_e\)
## @end tex
## @ifnottex
## @code{ye}
## @end ifnottex
## is not larger than
## @tex
## \(s_T \max(|y_f-y|)\).
## @end tex
## @ifnottex
## @code{sT*max(|yf-y|)}.
## @end ifnottex
## @item SettlingTime
## Time after which the absolute error
## @tex
## \(y_e\)
## @end tex
## @ifnottex
## @code{ye}
## @end ifnottex
## is not larger than
## @tex
## \(s_T |y_f-y_i|\).
## @end tex
## @ifnottex
## @code{sT*|yf-yi|}.
## @end ifnottex
## @item SettlingMin
## The minimum value of
## @tex
## \(y\)
## @end tex
## @ifnottex
## @code{y}
## @end ifnottex
## after it has risen (reached
## @tex
## \(r_U y_f\)).
## @end tex
## @ifnottex
## @code{rU*yf}).
## @end ifnottex
## @item SettlingMax
## The maximum value of
## @tex
## \(y\)
## @end tex
## @ifnottex
## @code{y}
## @end ifnottex
## after it has risen (reached
## @tex
## \(r_U y_f\)).
## @end tex
## @ifnottex
## @code{rU*yf}).
## @end ifnottex
## @item Overshoot
## Relative overshoot with respect to the normalized step response
## @tex
## \(y_n = (y - y_i)/(y_f - y_i)\).
## @end tex
## @ifnottex
## @code{yn = (y - yi)/(yf - yi)}.
## @end ifnottex
## @item Undershoot
## Relative undershoot with respect to the normalized step response
## @tex
## \(y_n = (y - y_i)/(y_f - y_i)\).
## @end tex
## @ifnottex
## @code{yn = (y - yi)/(yf - yi)}.
## @end ifnottex
## @item Peak
## Peak value of the absolute step response with respect to
## @tex
## \(y_i\).
## @end tex
## @ifnottex
## @code{yi}.
## @end ifnottex
## @item PeakTime
##  Time at which the peak value occurs.
## @end table
## @end table
##
## @seealso{step}
## @end deftypefn

function info = stepinfo (sys, varargin)

  ## Check input arguments

  if (nargin < 1)
    error ("stepinfo: at least one input argument required\n");
  endif

  if (! isa (sys, "lti"))
    error ("stepinfo: first argument must be an lti system\n");
  endif

  p = inputParser ();
  p.FunctionName = "stepinfo";
  vld_rtlimits = @(x) is_real_matrix (x) && ...
                      length (x(:)) == 2 && ...
                      all(x< 1) && all(x>0) && x(2) > x(1);
  p.addParameter ("RiseTimeLimits", [0.1 0.9],vld_rtlimits);
  vld_setthresh = @(x) isreal (x) && isscalar (x) && x < 1 && x > 0;
  p.addParameter ("SettlingTimeThreshold", 0.02,vld_setthresh);

  p.parse (varargin{:});

  ## Thresholds

  thresholds = struct ();
  thresholds.t_rise = p.Results.RiseTimeLimits;
  thresholds.t_settling = p.Results.SettlingTimeThreshold;
  thresholds.t_transient = p.Results.SettlingTimeThreshold;

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
    y_risen = cellfun (@(yf,yr) [thresholds.t_rise(2)*yf; yr], y_final, y_risen, "uniformoutput", false);
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

