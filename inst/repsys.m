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
## @deftypefn {Function File} {@var{rsys} =} repsys (@var{sys}, @var{m}, @var{n})
## @deftypefnx {Function File} {@var{rsys} =} repsys (@var{sys}, [@var{m}, @var{n}])
## @deftypefnx {Function File} {@var{rsys} =} repsys (@var{sys}, @var{m})
## Form a block transfer matrix of @var{sys} with @var{m} copies vertically
## and @var{n} copies horizontally.  If @var{n} is not specified, it is set to @var{m}.
## @code{repsys (sys, 2, 3)} is equivalent to @code{[sys, sys, sys; sys, sys, sys]}.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: May 2014
## Version: 0.1

function sys = repsys (varargin)

  if (nargin == 0)
    print_usage ();
  endif

  sys = repmat (varargin{:});   # repmat is overloaded for LTI systems

endfunction
