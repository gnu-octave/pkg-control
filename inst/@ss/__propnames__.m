## Copyright (C) 2009   Lukas F. Reichlin
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{propv}, @var{asgnvalv}] =} __propnames__ (@var{sys})
## @deftypefnx {Function File} {[@var{propv}, @var{asgnvalv}] =} __propnames__ (@var{sys}, @var{"specific"})
## Return the list of properties as well as the assignable values for a ss object sys.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.1

function [propv, asgnvalv] = __propnames__ (sys, flg)

  ## cell vector of ss-specific properties
  propv = {"a";
           "b";
           "c";
           "d";
           "stname"};

  ## cell vector of ss-specific assignable values
  asgnvalv = {"n-by-n matrix (n = number of states)";
              "n-by-m matrix (m = number of inputs)";
              "p-by-n matrix (p = number of outputs)";
              "p-by-m matrix";
              "n-by-1 cell vector of strings"};

  if (nargin == 1)
    [ltipropv, ltiasgnvalv] = __propnames__ (sys.lti);

    propv = [propv;
             ltipropv];

    asgnvalv = [asgnvalv;
                ltiasgnvalv];
  endif

endfunction