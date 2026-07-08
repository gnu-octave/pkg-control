## Copyright (C) 2026  Mitchell Thompkins <mitchell.thompkins@pm.me>
##
## This file is part of the control package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not,
## see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## Return the list of keys for a zpk object.

## Author: Mitchell Thompkins <mitchell.thompkins@pm.me>
## Created: June 2026
## Version: 0.1

function [keys, vals] = __sys_keys__ (sys, aliases = false)

  keys = {"z";
          "p";
          "k"};

  vals = {"p-by-m cell array of zero vectors (m = number of inputs)";
          "p-by-m cell array of pole vectors (p = number of outputs)";
          "p-by-m real-valued gain matrix"};

  if (aliases)
    ka = {"zeros"; "poles"; "gain"};
    keys = [keys; ka];
  endif

endfunction
