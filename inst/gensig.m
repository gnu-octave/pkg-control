## Copyright (C) 2009 Lukas Reichlin <lukas.reichlin@swissonline.ch>
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn{Function File} {[@var{u}, @var{t}] =} gensig (@var{sigtype}, @var{tau})
## @deftypefnx{Function File} {[@var{u}, @var{t}] =} gensig (@var{sigtype}, @var{tau}, @var{tfinal})
## @deftypefnx{Function File} {[@var{u}, @var{t}] =} gensig (@var{sigtype}, @var{tau}, @var{tfinal}, @var{tsam})
## Generate periodic signal. Useful in combination with lsim.
##
## @strong{Inputs}
## @table @var
## @item sigtype = "sin"
## Sine wave.
## @item sigtype = "cos"
## Cosine wave.
## @item sigtype = "square"
## Square wave.
## @item sigtype = "pulse"
## Periodic pulse.
## @item tau
## Duration of one period in seconds.
## @item tfinal
## Optional duration of the signal in seconds. Default duration is 5 periods
## @item tsam
## Optional sampling time in seconds. Default spacing is tau/64.
## @end table
##
## @strong{Outputs}
## @table @var
## @item u
## Vector of signal values.
## @item t
## Time vector of the signal.
## @end table
##
## @seealso{lsim}
## @end deftypefn

## Version: 0.1

function [u, t] = gensig (sigtype, tau, tfinal, tsam)
  
  if (nargin < 2 || nargin > 4)
    print_usage ();
  endif
  
  if (! ischar (sigtype))
    error ("gensig: first argument must be a string");
  endif
  
  if (nargin < 3)
    tfinal = 5 * tau;
  endif
  
  if (nargin < 4)
    tsam = tau / 64;
  endif
  
  t = (0 : tsam : tfinal)';
  
  sigtype = lower (sigtype);
  
  switch (sigtype(1:2))
    case "si"
      u = sin (2*pi/tau * t);
    case "co"
      u = cos (2*pi/tau * t);
    case "sq"
      u = +(rem (t, tau) >= tau/2);
    case "pu"
      u = +(rem (t, tau) < (1 - 1000*eps) * tsam);
    otherwise
      error ("gensig: invalid signal type");
  endswitch
  
endfunction

## FIXME: Add a test