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
## @deftypefn {Function File} {@var{bool} =} hasdelay (@var{sys})
## Return true if @var{sys} has a nonzero @var{InputDelay}, @var{OutputDelay},
## or @var{IODelay}.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} system.
## @end table
##
## @strong{Outputs}
## @table @var
## @item bool
## True if @var{sys} has any nonzero delay, false otherwise.
## @end table
## @end deftypefn

## Author: Prateek Ganguli <prateek.ganguli@gmail.com>
## Created: July 2026
## Version: 0.1

function bool = hasdelay (sys)

  if (nargin != 1 || ! isa (sys, "lti"))
    print_usage ();
  endif

  [indelay, outdelay, iodelay] = get (sys, "inputdelay", "outputdelay", "iodelay");

  bool = any (indelay(:) != 0) || any (outdelay(:) != 0) || any (iodelay(:) != 0);

endfunction


%!test
%! h = tf (10, [1 3 10]);
%! assert (hasdelay (h), false);

%!test
%! h = tf (10, [1 3 10], "IODelay", 0.25);
%! assert (hasdelay (h), true);

%!test
%! h = tf ([1 -1], [1 4 5], "InputDelay", 0.3);
%! assert (hasdelay (h), true);

%!error (hasdelay ())
%!error (hasdelay (1))
