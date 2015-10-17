## Copyright (C) 2009-2015   Lukas F. Reichlin
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
## Convert to descriptor state-space model.
## Since it is impossible to call a function inside @NAME without any
## argument of class NAME (except the constructor @NAME/NAME.m),
## the "constructor" for dss models is located outside folder @ss.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2010
## Version: 0.1

function sys = dss (sys)

  if (nargin != 1)
    print_usage ();
  endif

  if (! isa (sys, "ss"))
    sys = ss (sys);
  endif

  a = __sys_data__ (sys);

  sys = __set__ (sys, "e", eye (size (a)));

  ## FIXME: dss (1, 2, 3, 4, 5, tf()),  dss (1, ss(4)),  ...
  ##        Solve this by removing @lti/dss.m  and moving
  ##        all the functionality to the "constructor" dss.m
  ##        outside the '@lti' folder.

endfunction
