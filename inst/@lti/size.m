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
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with LTI Syncope.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {@var{nvec} =} size (@var{ltisys})
## @deftypefnx {Function File} {@var{n} =} size (@var{ltisys}, @var{idx})
## @deftypefnx {Function File} {[@var{ny}, @var{nu}] =} size (@var{ltisys})
## LTI model size, i.e. number of outputs and inputs.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.1

function [n, varargout] = size (ltisys, idx = 0)

  if (nargin > 2)
    print_usage ();
  endif

  ny = numel (ltisys.outname);  # WARNING: system matrices may change without
  nu = numel (ltisys.inname);   #          being noticed by the i/o names!

  switch (idx)
    case 0
      switch (nargout)
        case 0
          if (ny == 1)
            stry = "";
          else
            stry = "s";
          endif          
          if (nu == 1)
            stru = "";
          else
            stru = "s";
          endif
          disp (sprintf ("LTI model with %d output%s and %d input%s.", ny, stry, nu, stru));
        case 1
          n = [ny, nu];
        case 2
          n = ny;
          varargout{1} = nu;
        otherwise
          print_usage ();
      endswitch
    case 1
      n = ny;
    case 2
      n = nu;
    otherwise
      print_usage ();
  endswitch

endfunction