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
## @item N
## Number of poles of the approximations.
## @item M
## Number of zeros of the approximation, If omitted, the number of zeros is
## the same as the number @var{N} of poles.
## @end table
##
## More then one approximation can be requested by providing the additional
## order pairs @var{n2}, @var{m2}, ..., @var{nk}, @var{mk}.
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
## If no output argument is requested, the step response of the approximatons
## with orders ## @var{N1}, @var{M1} to @var{NL}, @var{ML} are plotted togehter
## with step delayed by the given dead-time @var{T}.
##
## When using the same numbers of poles and zeros (@var{m} = @var{n}), the step response of the
## resulting approximation shows a step at t = 0 which is untypical for a
## pure teim delay. Therefore, [1] proposes an approximaton with less zeros
## than poles (@var{m} < @var{n}).
##
## Algorithm based on:@*
## [1] Vajta, M. (2000). Some remarks on Pade-approximations, pp53-58,
##     Paper presented at 3rd TEMPUS-INTCOM Symposium on Intelligent Systems
##     in Control and Measurements 2000, Veszprém, Hungary.
## @end deftypefn

function varargout = pade (T, n, varargin)

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



  ## Parameter tests
  ## TODO

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
        sys = cell (1,N)
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

  t_end = 3*T;
  delay_text = sprintf ('Delay %10.3e s', T);

  t_delay = [0, T, T, t_end];
  y_delay = [0, 0, 1, 1];

  clf ();

  plot (t_delay, y_delay, 'linewidth', 1.5)

  hold on;

  leg = cell (1,N+1);
  leg{1} = delay_text;

  for i = 1:N
    [y,t] = step (tf (num{i}, den{i}), t_end);
    plot (t, y, 'linewidth', 1.5)
    leg{i+1} = ['Pade n=', num2str(length(den{i})-1), ', m=', num2str(length(num{i})-1)];
  endfor

  box on;
  grid on;

  legend (leg);
  legend ('location', 'southeast')

endfunction
