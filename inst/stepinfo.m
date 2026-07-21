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
## Array with two entries with relative values of
## @tex
## \(y_f\)
## @end tex
## @ifnottex
## @code{yf}
## @end ifnottex
## between which the rise time is determined.
## The values have to be between 0 and 1. The default is
## @tex
## \([r_L \, r_U] = [0.1 \, 0.9]\)
## @end tex
## @ifnottex
## @code{[rL rU] = [0.1 0.9]}
## @end ifnottex
## if this property is omitted.
## @item SettlingTimeThreshold
## Relative value of the absolute maximum or initial difference to
## @tex
## \(y_f\)
## @end tex
## @ifnottex
## @code{yf}
## @end ifnottex
## defining a tolerance range around
## @tex
## \(y_f\)
## @end tex
## @ifnottex
## @code{yf}
## @end ifnottex
## for TransientTime or SettlingTime (see below). The default is
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
## Time at which the peak value occurs.
## @item FinalValue
## Final value
## @tex
## \(y(t\to\infty)\)
## @end tex
## @ifnottex
## @code{y(t->Inf)}
## @end ifnottex
## of the step response, corresponds to the static gain.
## @item InitialJump
## Output value
## @tex
## \(y(t=0^+)\)
## @end tex
## @ifnottex
## @code{y(t=0+)}
## @end ifnottex
## immediately after the input step.
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

  p = inputParser ();   # Use Octvae's input parser for the input arguments
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

  ## Collect all system channels as separate systems in a cell array
  tfcell = cellfun (@(n,d) tf(n,d,sys.ts), num, den, "uniformoutput", false);

  ## Get stability (for each channel) and continuous-/discrete-time
  stability = cellfun (@isstable, tfcell);
  ct = isct (sys);

  ## Use a special version of __time_response__, with a longer time horizon
  sys_cell = cell ();
  sys_cell{1} = sys;    # __time_response__ expects a cell with systems
  [y,t] = __time_response__ ("stepinfo", sys_cell, 2);
  t = t{1,1};           # and returns a cell array
  y = y{1,1};
  N = length (t);

  ## Make a cell with all system outputs y
  y_cell = mat2cell (reshape (y,p*N,m), N*ones(1,p), ones (1,m));

  ## Get the static gain, y_final and y_init as cell arrays
  K = dcgain (sys);
  y_final = K;
  y_init = 0 * y_final; # for now, only y_init = 0 is considered
  y_final = mat2cell (y_final, ones(1,p), ones (1,m));
  y_init = mat2cell (y_init, ones(1,p), ones (1,m));

  ## Create structure array with all values

  info = struct ("RaiseTime", cell(p,m));

  for iy = 1:p
    for iu = 1:m

      ## Get values of current channel (from u(iu) to y(iy))
      y  = y_cell{iy,iu};

      if (! stability(iy,iu))

        ## Unstable system, output warning and set appropriate values

        warning ("stepinfo: system from u%d to y%d is not stable\n", iu, iy);

        rise_time = NaN;
        trans_time = NaN;
        settl_time = NaN;
        settl_min = NaN;
        settl_max = NaN;
        overshoot = NaN;
        undershoot = NaN;
        peak = Inf;
        peak_time = Inf;

      else

        ## Stable system

        ## Get values of current channel (from u(iu) to y(iy))
        yf = y_final{iy,iu};
        yi = y_init{iy,iu};

        ## Rise time
        t_rise = __get_rise_time__ (y, yf, yi, ct, t, thresholds.t_rise);
        rise_time = diff (t_rise);


        ## Transient time
        trans_time = __get_transient_time__ (y, yf, yi, ct, t, 't', thresholds.t_settling);

        ## Settling time
        settl_time = __get_transient_time__ (y, yf, yi, ct, t, 's', thresholds.t_settling);

        ## SettlingMin/Max
        risen_idx = find(t >= t_rise(2))(1);   # Time of upper risen level
        t_risen = t(risen_idx:end);            # Time vector from this time until end
        y_risen = y(risen_idx:end);            # Output from this time until end

        if (ct)
          y_risen = [thresholds.t_rise(2)*yf; y_risen];
        endif

        settl_min = min (y_risen);
        settl_max = max (y_risen);

        ## Overshoot / Undershoot
        y_rel = y - yi;
        y_norm = y_rel./(yf - yi);

        overshoot = max (0, 100 * max (y_norm-1));
        undershoot = max (0, -100 * min (y_norm));

        ## Peak and Peak Time
        s = sign (y(end) - yi); # don't use y_final, which is 0 for D systems
        if (s == 0)
          ## In case of D-systems, y_final does not
          s = y(1) - yi;
        endif
        y_peak = s * y_rel;

        [peak idx] = max (y_peak);
        peak_time = t(idx);

      endif

    info = setfield (info, {iy,iu}, "RiseTime" , rise_time);
    info = setfield (info, {iy,iu}, "TransientTime" , trans_time);
    info = setfield (info, {iy,iu}, "SettlingTime" , settl_time);
    info = setfield (info, {iy,iu}, "SettlingMin" , settl_min);
    info = setfield (info, {iy,iu}, "SettlingMax" , settl_max);
    info = setfield (info, {iy,iu}, "Overshoot" , overshoot);
    info = setfield (info, {iy,iu}, "Undershoot" , undershoot);
    info = setfield (info, {iy,iu}, "Peak" , peak);
    info = setfield (info, {iy,iu}, "PeakTime" , peak_time);
    info = setfield (info, {iy,iu}, "StaticGain" , yf);
    info = setfield (info, {iy,iu}, "InitialJump" , y(1));

    endfor
  endfor


endfunction



function t_rise = __get_rise_time__ (y, y_final, y_init, ct, t, t_tol)

  y_rise = arrayfun (@(x) t_tol * x, y_final, "uniformoutput", false);
  t_rise = [-1 -2];

  s = sign (y_final - y_init);
  if (s == 0)
    s = sign (y(1) - y_init);
  endif

  for i = 1:2

    j = (find (s*y > s*y_rise{1}(i)))(1);

    if (! ct)
      t_rise(i) = t(j);
      continue
    endif

    if (j > 1)
      t_rise(i) = interp1 ([y(j-1) y(j)], [t(j-1) t(j)], y_rise{1}(i));
    else
      t_rise(i) = 0;
    endif

  endfor

endfunction




function t_transient = __get_transient_time__ (y, y_final, y_init, ct, t, t_type, t_tol)

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

