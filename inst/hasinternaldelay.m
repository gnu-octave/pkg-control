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
## @deftypefn {Function File} {@var{bool} =} hasinternaldelay (@var{sys})
## Return true if @var{sys} has a nonzero @var{InternalDelay}.
## @var{InternalDelay} is a @code{ss}-specific property; any other
## @acronym{LTI} object (@code{tf}, @code{zpk}) always returns false.
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
## True if @var{sys} has any nonzero internal delay, false otherwise.
## @end table
## @end deftypefn

## Author: Prateek Ganguli <prateek.ganguli@gmail.com>
## Created: July 2026
## Version: 0.1

function bool = hasinternaldelay (sys)

  if (nargin != 1 || ! isa (sys, "lti"))
    print_usage ();
  endif

  if (! isa (sys, "ss"))
    bool = false;
    return;
  endif

  internaldelay = get (sys, "internaldelay");

  bool = any (internaldelay(:) != 0);

endfunction


%!test
%! h = tf (10, [1 3 10]);
%! assert (hasinternaldelay (h), false);

%!test
%! h = zpk ([], -2, 3);
%! assert (hasinternaldelay (h), false);

%!test
%! h = ss (-1, 1, 1, 0);
%! assert (hasinternaldelay (h), false);

%!test
%! h = ss (-1, 1, 1, 0);
%! h = set (h, "internaldelay", 0.5);
%! assert (hasinternaldelay (h), true);

%!error (hasinternaldelay ())
%!error (hasinternaldelay (1))
