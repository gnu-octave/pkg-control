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
## Class constructor.  For internal use only.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.2

function p = tfpoly (a)

  superiorto ("double");

  switch (nargin)
    case 0
      p = struct ("poly", []);
      p = class (p, "tfpoly");

    case 1
      if (isa (a, "tfpoly"))
        p = a;
        return;
      elseif (is_vector (a))
        if (! is_real_vector (a))
          warning (["The resulting transfer function might have complex coefficients.\n", ...
           "         Most functions from the Control Package won't work for transfer\n", ...
           "         functions with complex coefficients\n"]);
        endif
        p.poly = reshape (a, 1, []);
        p = class (p, "tfpoly");
        p = __remove_leading_zeros__ (p);
      else
        error ("tfpoly: argument must be a vector");
      endif

    otherwise
      print_usage ();

  endswitch

endfunction
