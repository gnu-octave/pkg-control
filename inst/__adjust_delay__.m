## Copyright (C) 2026        Prateek Ganguli
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
## Check whether a delay value has the required size, expanding a scalar
## to the full size if necessary, and enforce non-negative, finite values
## (integer-valued for discrete-time models).  Used by lti/set and __set__.

## Author: Prateek Ganguli <prateek.ganguli@gmail.com>
## Created: July 2026
## Version: 0.1

function val = __adjust_delay__ (val, req_size, isdiscrete, caller = "lti: set")

  if (! is_real_matrix (val))
    error ("%s: delay value must be a real-valued matrix", caller);
  endif

  if (isscalar (val))
    val = val * ones (req_size);
  elseif (! isequal (size (val), req_size))
    error ("%s: delay matrix must have size %dx%d", caller, req_size(1), req_size(2));
  endif

  if (any (val(:) < 0) || ! all (isfinite (val(:))))
    error ("%s: delay values must be non-negative and finite", caller);
  endif

  if (isdiscrete && any (val(:) != fix (val(:))))
    error ("%s: delay values must be integer-valued for discrete-time models", caller);
  endif

endfunction


%!assert (__adjust_delay__ (0.3, [1, 1], false), 0.3)
%!assert (__adjust_delay__ (0.3, [2, 1], false), [0.3; 0.3])
%!assert (__adjust_delay__ ([0.1; 0.2], [2, 1], false), [0.1; 0.2])
%!assert (__adjust_delay__ (0.25, [2, 2], false), [0.25, 0.25; 0.25, 0.25])
%!assert (__adjust_delay__ (3, [1, 1], true), 3)
%!assert (__adjust_delay__ (0, [1, 1], true), 0)

%!error (__adjust_delay__ (0.3, [1, 1], true))
%!error (__adjust_delay__ (-0.1, [1, 1], false))
%!error (__adjust_delay__ ([0.1, 0.2], [2, 1], false))
%!error (__adjust_delay__ ("a", [1, 1], false))
%!error (__adjust_delay__ (Inf, [1, 1], false))
