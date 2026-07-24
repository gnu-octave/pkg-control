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
## @deftypefn {Function File} {@var{delay} =} totaldelay (@var{sys})
## Return the p-by-m matrix of total delay (in seconds for continuous-time
## models, samples for discrete-time models) from each input to each
## output, i.e. @code{IODelay + InputDelay + OutputDelay}.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} system.
## @end table
##
## @strong{Outputs}
## @table @var
## @item delay
## p-by-m real matrix.  @code{delay(i,j)} is the total delay from
## input j to output i.
## @end table
## @end deftypefn

## Author: Prateek Ganguli <prateek.ganguli@gmail.com>
## Created: July 2026
## Version: 0.1

function delay = totaldelay (sys)

  if (nargin != 1 || ! isa (sys, "lti"))
    print_usage ();
  endif

  [indelay, outdelay, iodelay] = get (sys, "inputdelay", "outputdelay", "iodelay");
  [p, m] = size (sys);

  delay = iodelay + repmat (outdelay, 1, m) + repmat (indelay.', p, 1);

endfunction


%!test
%! h = tf (10, [1 3 10], "IODelay", 0.25);
%! assert (totaldelay (h), 0.25);

%!test
%! h = tf ({1, 1; 1, 1}, {[1 1], [1 2]; [1 3], [1 4]});
%! h.InputDelay = [0.1; 0.2];
%! h.OutputDelay = [0.3; 0.4];
%! assert (totaldelay (h), [0.4, 0.5; 0.5, 0.6], 1e-10);

%!error (totaldelay ())
%!error (totaldelay (1))
