## Copyright (C) 2023-2024 Torsten Lilge
##
## This file is part of the Control package for GNU Octave
##
## This is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This software is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with the Control package. If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {@var{st} =} doc_control (@var{fcn1}, @var{fcn2}, ...)
## Open online documentation of the Control package in the system's
## standard browser.
##
## @strong{Inputs}
## @table @var
## @item fcn1, fcn2, ...
## Function names for which the documentation should be displayed.
## If no function name is given, the index of the documentation is shown.
## If one of the function names is 'license', copyright and license
## information are displayed.
## @end table
##
## @strong{Outputs}
## @table @var
## @item st
## The return value @var{st} has one of the values:
## @itemize @bullet
## @item
## @samp{0} Found and opened system browser successfully.
## @item
## @samp{1} System browser not found.
## @item
## @samp{2} System browser found, but an error occurred.
## @end itemize
## @end table
##
## @end deftypefn

## Author: Torsten Lilge <ttl-octave@mailbox.org>
## Created: December 2023
## Version: 0.1

function st = doc_control (varargin)

  base = 'https://gnu-octave.github.io/pkg-control/';
  license = 'https://github.com/gnu-octave/pkg-control/blob/main/COPYING';

  status = 0;

  if nargin == 0

    status = web (base);

  else

    for i = 1:nargin

      if ischar (varargin{i})
        fcn = strtrim (varargin{i});
        if (strcmp (fcn, 'license'))
          url = license;
        else
          url = [base, strrep(fcn, '/', '_'), '.html'];
        endif
        sti = web (url);
        if sti == 1
          status = sti;
          break;
        elseif status == 0
          status = sti;
        endif
      else
        warning ("argument %d is not a string\n", i);
      endif

    endfor

    if nargout > 0
      st = status;
    endif

  endif

endfunction
