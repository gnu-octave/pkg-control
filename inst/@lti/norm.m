## Copyright (C) 2009   Lukas F. Reichlin
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
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {@var{gain} =} norm (@var{sys}, @var{2})
## @deftypefnx {Function File} {[@var{gain}, @var{wpeak}] =} norm (@var{sys}, @var{inf})
## @deftypefnx {Function File} {[@var{gain}, @var{wpeak}] =} norm (@var{sys}, @var{inf}, @var{tol})
## Return H-2 or L-inf norm of LTI model.
## Uses SLICOT AB13BD and AB13DD by courtesy of NICONET e.V.
## <http://www.slicot.org>
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: November 2009
## Version: 0.2

function [gain, varargout] = norm (sys, ntype = "2", tol = 0.01)

  if (nargin > 3)  # norm () is caught by built-in function
    print_usage ();
  endif

  if (isnumeric (ntype) && isscalar (ntype))
    if (ntype == 2)
      ntype = "2";
    elseif (isinf (ntype))
      ntype = "inf";
    else
      error ("lti: norm: invalid norm type");
    endif
  elseif (ischar (ntype))
    ntype = lower (ntype);
  else
    error ("lti: norm: invalid norm type");
  endif

  switch (ntype)
    case "2"
      gain = h2norm (sys);

    case "inf"
      [gain, varargout{1}] = linfnorm (sys, tol);

    otherwise
      error ("lti: norm: invalid norm type");
  endswitch

endfunction


function gain = h2norm (sys)

  if (isstable (sys))
    [a, b, c, d] = ssdata (sys);
    digital = ! isct (sys);
    
    if (! digital && ! all (all (d == 0)))  # continuous and non-zero feedthrough
      gain = inf;
    else
      [gain, iwarn] = slab13bd (a, b, c, d, digital);
    
      if (iwarn)
        warning ("lti: norm: slab13bd: iwarn = %d", iwarn);
      endif
    endif
  else
    gain = inf;
  endif

endfunction


function [gain, wpeak] = linfnorm (sys, tol = 0.01)

  [a, b, c, d, tsam] = ssdata (sys);
  digital = ! isct (sys);
  tol = max (tol, 100*eps);
  
  [fpeak, gpeak] = slab13dd (a, b, c, d, digital, tol);
  
  if (fpeak(2) > 0)
    if (digital)
      wpeak = fpeak(1) / tsam;
    else
      wpeak = fpeak(1);
    endif
  else
    wpeak = inf;
  endif
  
  if (gpeak(2) > 0)
    gain = gpeak(1);
  else
    gain = inf;
  endif

endfunction
