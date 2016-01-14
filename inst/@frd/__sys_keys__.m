## Copyright (C) 2009-2016   Lukas F. Reichlin
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
## @deftypefn {Function File} {[@var{keys}, @var{vals}] =} __sys_keys__ (@var{sys})
## @deftypefnx {Function File} {[@var{keys}, @var{vals}] =} __sys_keys__ (@var{sys}, @var{aliases})
## Return the list of keys as well as the assignable values for a frd object sys.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2010
## Version: 0.3

function [keys, vals] = __sys_keys__ (sys, aliases = false)

  ## cell vector of frd-specific keys
  keys = {"H";
          "w"};

  ## cell vector of frd-specific assignable values
  vals = {"p-by-m-by-l array of complex frequency responses";
          "l-by-1 vector of real frequencies (l = length (w))"};

  if (aliases)
    ka = {"response";
          "frequency"};
    keys = [keys; ka];
  endif

endfunction
