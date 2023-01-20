## Copyright (C) 2022  Torsten Lilge
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
## @deftypefn {Function File} {@var{bool} =} isstaticgain (@var{sys})
## Determine whether @acronym{LTI} model is a static gain.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} system.
## @end table
##
## @strong{Outputs}
## @table @var
## @item bool = 0
## @var{sys} is a dynamical system
## @item bool = 1
## @var{sys} is a static gain
## @end table
## @end deftypefn

## Author: Torsten Lilge <ttl-octave@mailbox.org>
## Created: October 2022
## Version: 0.1

function static_gain = isstaticgain (ltisys)

  if (nargin == 0)
    print_usage ();
  endif

  num = ltisys.num;
  den = ltisys.den;

  static_gain = false;

  if (isempty (den))                    # tf (num, []),  where [] could be {} as well
    if (isempty (num))                  # tf ([], [])
      num = den = {};
      static_gain = true;
    elseif (is_real_matrix (num))       # static gain  tf (matrix),  tf (matrix, [])
      num = num2cell (num);
      den = num2cell (ones (size (num)));
      static_gain = true;
    endif
  endif

  if (! iscell (num))
    num = {num};
  endif
  if (! iscell (den))
    den = {den};
  endif

  ## Now check for static gain if all tfs have size num and size den of one
  num_scalar = cellfun (@(p) (find (p != 0, 1) == length (p)) || (length (find (p != 0, 1)) == 0), num);
  den_scalar = cellfun (@(p) (find (p != 0, 1) == length (p)) || (length (find (p != 0, 1)) == 0), den);
  if (all (num_scalar) && all (den_scalar))
    ## All tf components are of the form b0/a0 (static gain)
    static_gain = true;
  endif

endfunction
