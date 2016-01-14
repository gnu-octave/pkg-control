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
## @deftypefn {Function File} {@var{dat} =} ifft (@var{dat})
## Transform iddata objects from frequency to time domain.
##
## @strong{Inputs}
## @table @var
## @item dat
## iddata set containing signals in frequency domain.
## The frequency values must be distributed equally from 0
## to the Nyquist frequency.  The Nyquist frequency is
## only included for even signal lengths.
## @end table
##
## @strong{Outputs}
## @table @var
## @item dat
## iddata identification dataset in time domain.
## In order to preserve signal power and noise level,
## the FFTs are normalized by multiplying each transform
## by the square root of the signal length.
## @end table
##
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: April 2012
## Version: 0.2

function dat = ifft (dat)

  if (nargin > 1)       # no need to test nargin == 0, this is handled by built-in ifft
    print_usage ();
  endif

  if (dat.timedomain)   # dat is always an iddata set, otherwise built-in ifft would be called
    return;
  endif

  if (any (cellfun (@(w) w(1) >= eps, dat.w)))
    error ("iddata: ifft: first frequency must be zero");
  endif

  if (any (cellfun (@(w) any (abs (diff (w, 2)) > 1e-4*w(2:end-1)), dat.w)))
    error ("iddata: ifft: require linearly spaced frequency vectors");
  endif

  [x, ~, ~, e] = size (dat);

  x = x(:);
  n = num2cell (x);
  nconj = num2cell (x - rem (x, 2));

  dat.y = cellfun (@(y, n, nconj) real (ifft ([y; conj(y(nconj:-1:2, :))], [], 1)) * sqrt (n+nconj), dat.y, n, nconj, "uniformoutput", false);
  dat.u = cellfun (@(u, n, nconj) real (ifft ([u; conj(u(nconj:-1:2, :))], [], 1)) * sqrt (n+nconj), dat.u, n, nconj, "uniformoutput", false);

  ## ifft (x, n, dim=1) because x could be a row vector (n=1)
  
  dat.w = {};   % dat.w = repmat ({[]}, e, 1); ???
  dat.timedomain = true;

endfunction
