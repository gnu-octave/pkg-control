## Copyright (C) 2011   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{nvec} =} size (@var{sys})
## @deftypefnx {Function File} {@var{n} =} size (@var{sys}, @var{dim})
## @deftypefnx {Function File} {[@var{n}, @var{p}, @var{m}, @var{ne}] =} size (@var{sys})
## LTI model size, i.e. number of outputs and inputs.
##
## @strong{Inputs}
## @table @var
## @item sys
## LTI system.
## @item dim
## If given a second argument, @command{size} will return the size of the
## corresponding dimension.
## @end table
##
## @strong{Outputs}
## @table @var
## @item nvec
## Row vector.  The first element is the number of outputs (rows) and the second
## element the number of inputs (columns).
## @item n
## Scalar value.  The size of the dimension @var{dim}.
## @item p
## Number of outputs.
## @item m
## Number of inputs.
## @end table
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2011
## Version: 0.1

function [ret, varargout] = size (dat, dim = 0)

  if (nargin > 2)
    print_usage ();
  endif

  p = numel (dat.outname);
  m = numel (dat.inname);
  n = rows (dat.y);

  switch (dim)
    case 0
      switch (nargout)
        case 0
          if (p == 1)
            stry = "";
          else
            stry = "s";
          endif          
          if (m == 1)
            stru = "";
          else
            stru = "s";
          endif
          disp (sprintf ("IDDATA set with %d output%s and %d input%s.", p, stry, m, stru));
        case 1
          ret = [n, p, m, 1];
        case {2, 3, 4}
          ret = n;
          varargout{1} = p;
          varargout{2} = m;
          varargout{3} = 1;
        otherwise
          print_usage ();
      endswitch
    case 1
      ret = n;
    case 2
      ret = p;
    case 3
      ret = m;
    case 4
      ret = 1;
    otherwise
      print_usage ();
  endswitch

endfunction
