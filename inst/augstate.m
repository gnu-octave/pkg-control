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
## @deftypefn {Function File} {@var{augsys} =} augstate (@var{sys})
## Append state vector x of system @var{sys} to output vector y.
##
## @example
## @group
## .                  .
## x = A x + B u      x = A x + B u
## y = C x + D u  =>  y = C x + D u
##                    x = I x + O u
## @end group
## @end example
##
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2015
## Version: 0.1

function augsys = augstate (sys)

  if (nargin != 1 || ! isa (sys, "lti"))
    print_usage ();
  endif

  if (! isa (sys, "ss"))
    warning ("augstate: system not in state-space form");
    sys = ss (sys);
  endif

  [a, b, c, d, e, tsam] = dssdata (sys, []);
  [inn, stn, outn, ing, outg] = get (sys, "inname", "stname", "outname", "ingroup", "outgroup");

  [p, m] = size (d);
  n = rows (a);

  caug = vertcat (c, eye (n));
  daug = vertcat (d, zeros (n, m));

  outname = vertcat (outn, stn);  
  
  augsys = dss (a, b, caug, daug, e, tsam);
  augsys = set (augsys, "inname", inn, "stname", stn, "outname", outname, ...
                        "ingroup", ing, "outgroup", outg);

endfunction
