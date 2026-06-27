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

  ## Get information on the system

  [p, m] = size (sys);
  [num, den] = tfdata (sys);
  tfcell = cellfun (@tf, num, den, "uniformoutput", false);

  stability = mat2cell (cellfun (@isstable, tfcell), ...
                        ones (1,p),  ones (1,m));

  K = dcgain (sys);
  yfinal = K;

  [y,t] = step (sys);
  N = length (t);

  y_cell = mat2cell (reshape (y,p*N,m), N*ones(1,p), ones (1,m));

  ## Prepare array of structure as output argument
  y_rise = arrayfun (@(x) [0.1 0.9] * x, yfinal, "uniformoutput", false);

  info = struct ("RiseTime", cell (size (stability)));


  t_rise = cellfun (@(yc,yr,st) __get_rise_time__ (yc,yr,st,t), ...
                    y_cell, y_rise, stability, "uniformoutput", false);
  info = struct ("RiseTime", ...
                 cellfun (@diff, t_rise, "uniformoutput", false));

endfunction



function t_rise = __get_rise_time__ (y, y_rise, stability, t)

  t_rise = [-1 -2];

  if (! stability)
    return;
  endif

  for i = 1:2
    j = (find (y > y_rise(i)))(1);
    if (j > 1)
      t_rise(i) = interp1 ([y(j-1) y(j)], [t(j-1) t(j)], y_rise(i));
    else
      if (i == 1)
        t_rise(i) = 0;
      else
        t_rise(i) = Inf;
      endif
    endif
  endfor

endfunction
