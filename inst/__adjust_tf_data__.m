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
## Common code for adjusting TF model data.
## Used by tf and __set__.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2015
## Version: 0.1

function [num, den, tsam, tfvar] = __adjust_tf_data__ (num, den, tsam = -1)

  static_gain = false;
  if (isempty (den))                    # tf (num, []),  where [] could be {} as well
    if (isempty (num))                  # tf ([], [])
      num = den = {};
      static_gain = true;
    elseif (is_matrix (num))       # static gain  tf (matrix),  tf (matrix, [])
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

  ## NOTE: the 'tfpoly' constructor checks its vector as well,
  ##       but its error message would make little sense for users
  ##       and would make it hard for them to identify the invalid argument.

  ## NOTE: this code is nicer, but there would be conflicts in  @tf/__set__.m
  ##       e.g.  sys = tf ([5, 22], [1, 2]),  sys.num = 5
  ##
  ## if (! is_real_vector (num{:}, 1))     # dummy argument 1 needed if num is empty cell
  ##   error ("tf: first argument 'num' requires a real-valued, non-empty vector or a cell of such vectors");
  ## endif
  ## if (! is_real_vector (den{:}, 1))
  ##  error ("tf: second argument 'den' requires a real-valued, non-empty vector or a cell of such vectors");
  ## endif
  ##
  ## num = cellfun (@tfpoly, num, "uniformoutput", false);
  ## den = cellfun (@tfpoly, den, "uniformoutput", false);

  try
    num = cellfun (@tfpoly, num, "uniformoutput", false);
  catch
    error ("tf: numerator 'num' must be a non-empty vector or a cell of such vectors");
  end_try_catch

  try
    den = cellfun (@tfpoly, den, "uniformoutput", false);
  catch
    error ("tf: denominator 'den' must be a real-valued, non-empty vector or a cell of such vectors");
  end_try_catch

  if (any (cellfun (@is_zero, den)(:)))
    error ("tf: denominator(s) cannot be zero");
  endif

  if (static_gain)
    tfvar = "x";
  elseif (tsam == 0)
    tfvar = "s";
  else
    tfvar = "z";
  endif

endfunction
