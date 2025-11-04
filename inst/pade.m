## Copyright (C) 2025 Torten Lilge
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
## @deftypefn {Function File} pade (@var{T}, @var{n})
## @deftypefnx {Function File} pade (@var{T}, @var{n}, @var{m})
## @deftypefnx {Function File} pade (@var{T}, @var{n1}, @var{m1}, @var{n2}, @var{m2}, ..., @var{nk}, @var{mk})
## @deftypefnx {Function File} @var{sys} = pade (...)
## @deftypefnx {Function File} @var{sys} = pade (...)
## @deftypefnx {Function File} [@var{num}, @var{den}] = pade (...)
## @deftypefnx {Function File} [@var{num}, @var{den}] = pade (...)
##
## Calculate Padé approximation of a dead-time by zeros and poles
##
## @strong{Inputs}
## @table @var
## @item T
## Dead-time to be approximated.
## @item n
## Number of poles of the approximations.
## @item m
## Number of zeros of the approximation, If omitted, the number of zeros is
## the same as the number @var{n} of poles.
## @end table
##
## More then one approximation can be requested by providing the pairs
## pairs @var{n1}, @var{m1}, @var{n2}, @var{m2}, ..., @var{nk}, @var{mk}.
##
## @strong{Outputs}
## @table @var
## @item sys
## LTI system with the poles and zeros of the Padé approximation. @var{sys}
## is given as transfer function. If more than one approximation is requested,
## @var{sys} is a cell array of LTI systems.
## @item num
## Numerator polynomial of the resulting transfer function.
## If more than one approximation is requested,
## @var{num} is a cell array of numerator polynomials.
## @item den
## Denominator polynomial of the resulting transfer function.
## If more than one approximation is requested,
## @var{den} is a cell array of denominator polynomials.
## @end table
##
## If no output argument is requested, the step response and the bode diagram
## of the approximatons with orders @var{n1}, @var{m1} to @var{nk}, @var{nk}
## are plotted togehter with step delayed by the given dead-time @var{T}.
##
## When using the same numbers of poles and zeros (@var{m} = @var{n}), the step response of the
## resulting approximation shows a step at t = 0 which is untypical for a
## pure time delay. However, it has a magnitude of 1 (0 dB) over all
## frequencies. In [1], an approximaton having less zeros than poles
## (@var{m} < @var{n}) is suggested as an alternative approach, resulting
## in a better step response but decresing magnitude for higher frequencies.
##
## Algorithm based on:@*
## [1] Vajta, M. (2000). Some remarks on Pade-approximations, pp53-58,
##     Paper presented at 3rd TEMPUS-INTCOM Symposium on Intelligent Systems
##     in Control and Measurements 2000, Veszprém, Hungary.
##
## @seealso{@@tf/tf}
## @end deftypefn

function varargout = pade (T, n, varargin)

  ## Parameter check
  if nargin < 2
    print_usage ();
  endif

  if (nargin > 2) && (mod (nargin,2) == 0)
    msg = ["pade: for more than one approximation, the number of input ",...
           "arguments must be odd (number of poles and zeros in pairs)\n"];
    error (msg);
  endif

  if ! (isreal (T) && T >= 0)
    error ("pade: dead time T must be real and non-negative\n");
  endif

  if (! __isinteger__ (n)) || ...
     (nargin > 2 && ! all (cellfun (@__isinteger__, varargin)))
    error ("pade: number of poles or zeros must be integer numbers\n");
  endif

  ## Get parameters
  if nargin > 2
    N = (nargin-1)/2;
    ni = zeros(1,(nargin-1)/2);
    mi = zeros(1,(nargin-1)/2);
    idx = 0;
    for i = 1:N
      if i == 1
        ni(i) = n;
      else
        ni(i) = varargin{++idx};
      endif
      mi(i) = varargin{++idx};
    endfor
  else
    N = 1;
    ni = n;
    mi = n;
  endif

  ## Check for causality
  if any (ni < mi)
    error ("pade: non-causal approximation requested (n < m)\n");
  endif


  ## Compute numeratr and denominator polinomials following [1], Section 3
  # T as array of correct length
  num = cell (1,N);
  den = cell (1,N);

  for i = 1:N

    n = ni(i);
    m = mi(i);

    Tm = T * ones (1,m+1);
    Tn = T * ones (1,n+1);
    # Powers in the polynomials
    pm = m:-1:0;
    pn = n:-1:0;
    # Factorials (constant 'c' and variable 'v' depending on power)
    fmc = factorial (m) / factorial (m+n);
    fmv = factorial (m+n-pm) ./ factorial (pm) ./ factorial (m-pm);
    fnc = factorial (m) / factorial (m+n);
    fnv = factorial (m+n-pn) ./ factorial (pn) ./ factorial (n-pn);

    # Numerator and denominator
    num{i} = fmc * (-Tm) .^ pm .* fmv;
    den{i} = fnc * Tn .^ pn .* fnv;
    # Adjust static gain to 1
    gain = num{i}(end)/den{i}(end);
    den{i} = den{i}*gain;

  endfor

  ## Determine what to return
  if nargout == 0

    __pade_plot__ (num, den, T);

  else

    if N == 1
      num = num{1};
      den = den{1};
    endif

    if nargout == 1

      if N == 1
        sys = tf (num, den);
      else
        sys = cell (1,N);
        for i = 1:N
          sys{i} = tf (num{i}, den{i});
        endfor
      endif

      varargout{1} = sys;

    endif

    if nargout == 2

      varargout{1} = num;
      varargout{2} = den;

    endif

  endif

endfunction


## Plot the delay and its Pade approximation
function __pade_plot__ (num, den, T)

  N = length(num);

  t_end = 2.5*T;
  if t_end == 0
    t_end = 1;
  endif

  delay_text = sprintf ('Delay %10.3e s', T);

  ## Step response
  t_delay = [0, T, T, t_end];
  y_delay = [0, 0, 1, 1];

  step_response_title = 'Pade approximation: Step response';
  bode_plot_title = 'Pade approximation: Bode plot';

  fig_handles = findobj ('Type', 'figure');

  h_step = __handle_from_name__ (fig_handles, step_response_title);
  h_bode = __handle_from_name__ (fig_handles, bode_plot_title);

  figure (h_step);
  clf ();
  set (gcf (), 'numbertitle', 'off');
  set (gcf (), 'name', step_response_title);

  hold on;

  sys = cell(1,N);
  leg = cell (1,N+1);
  leg{1} = delay_text;

  for i = 1:N
    sys{i} = tf (num{i}, den{i});
    [y,t] = step (sys{i}, t_end);
    plot (t, y, 'linewidth', 1.5)
    leg{i} =  ['Pade n=', num2str(length(den{i})-1), ', m=', num2str(length(num{i})-1)];
  endfor

  plot (t_delay, y_delay, '-.k', 'linewidth', 1)
  leg{end} = delay_text;

  box on;
  grid on;

  legend (leg);
  legend ('location', 'southeast')

  ## Bode plot
  figure (h_bode);
  clf ();
  set (gcf (), 'numbertitle', 'off');
  set (gcf (), 'name', bode_plot_title);
  bode (sys{:});

  xlim10 = log10 (xlim ());

  w = logspace (xlim10(1),xlim10(2),1000);
  pha_delay = -w*T/pi*180;

  hc = get (h_bode, 'children');
  amp = hc(2);
  pha = hc(4);

  hold (amp, 'on');
  plot (amp, w, pha_delay, '-.k', 'linewidth', 1);
  hold (amp, 'off');
  legend ('location', 'southwest');

  hold (pha, 'on');
  plot (pha, w, zeros (size(w)), '-.k', 'linewidth', 1);  # zeros, because dB
  hold (pha, 'off');

  amp_leg = get(amp, "__legend_handle__");
  set (amp_leg, "string", leg);
  set (amp_leg, 'location', 'southwest');

  pha_leg = get(pha, "__legend_handle__");
  set (pha_leg, "string", leg);
  set (pha_leg, 'location', 'southwest');

endfunction


## Get figure handles from name
function h = __handle_from_name__ (handles, name);

  idx = find (strcmp (get (handles, 'name'), name));
  if (isempty (idx))
    h = figure ();  # not found, create new figure
  else
    h = handles (idx(1));   # name found, get corresponding handle
  endif

endfunction


## Check for integer
function intx = __isinteger__ (x)

  intx = round (x) == x;

endfunction


%!demo
%! pade (1,4,4,4,3,4,2)

## Examples from [1, Table 1]

%!test
%! T = 2*rand ();
%! [num33,den33] = pade (T,3);
%! num33e = [ -T^3   12*T^2  -60*T  120 ];
%! den33e = [  T^3   12*T^2   60*T  120 ];
%! num33 = num33/num33(end);
%! den33 = den33/num33(end);
%! den33e = den33e/num33e(end);
%! num33e = num33e/num33e(end);
%! assert (num33, num33e, 10-4);
%! assert (den33, den33e, 10-4);

%!test
%! T = 3*rand ();
%! [num43,den43] = pade (T,4,3);
%! num43e = [         -T^3   60*T^2 -360*T  840 ];
%! den43e = [  T^4  16*T^3  120*T^2  480*T  840 ];
%! den43 = den43/num43(end);
%! num43 = num43/num43(end);
%! den43e = den43e/num43e(end);
%! num43e = num43e/num43e(end);
%! assert (num43, num43e, 10-4);
%! assert (den43, den43e, 10-4);

