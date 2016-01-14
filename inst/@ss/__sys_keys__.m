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
## @deftypefn {Function File} {[@var{keys}, @var{vals}] =} __sys_keys__ (@var{sys})
## @deftypefnx {Function File} {[@var{keys}, @var{vals}] =} __sys_keys__ (@var{sys}, @var{aliases})
## Return the list of keys as well as the assignable values for a ss object sys.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: September 2009
## Version: 0.3

function [keys, vals] = __sys_keys__ (sys, aliases = false)

  ## cell vector of ss-specific keys
  keys = {"a";
          "b";
          "c";
          "d";
          "e";
          "stname";
          "scaled"};

  ## cell vector of ss-specific assignable values
  vals = {"n-by-n matrix (n = number of states)";
          "n-by-m matrix (m = number of inputs)";
          "p-by-n matrix (p = number of outputs)";
          "p-by-m matrix";
          "n-by-n matrix";
          "n-by-1 cell vector of strings";
          "scalar logical value"};

  if (aliases)
    ka = {"statename"};
    keys = [keys; ka];
  endif

endfunction
