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
## @deftypefn {Function File} {[@var{keys}, @var{vals}] =} __iddata_keys__ (@var{dat})
## Return the list of keys as well as the assignable values for an iddata set.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: February 2012
## Version: 0.2

function [keys, vals] = __iddata_keys__ (dat, aliases = false)

  [n, p, m, e] = size (dat);

  ## cell vector of iddata-specific keys
  keys = {"y";
          "outname";
          "outunit";
          "u";
          "inname";
          "inunit";
          "tsam";
          "timeunit";
          "expname";
          "name";
          "notes";
          "userdata"};

  ## cell vector of iddata-specific assignable values
  vals = {sprintf("(%dx1) cell vector of (nx%d) matrices", e, p);
          sprintf("(%dx1) cell vector of strings", p);
          sprintf("(%dx1) cell vector of strings", p);
          sprintf("(%dx1) cell vector of (nx%d) matrices", e, m);
          sprintf("(%dx1) cell vector of strings", m);
          sprintf("(%dx1) cell vector of strings", m);
          sprintf("(%dx1) cell vector of scalars", e);
          "string";
          sprintf("(%dx1) cell vector of strings", e);
          "string";
          "string or cell of strings";
          "any data type"};

  if (aliases)
    ka = {"outdata";
          "outputdata";
          "outputname";
          "outputunit";
          "indata";
          "inputdata";
          "inputname";
          "inputunit";
          "experimentname";
          "w";
          "frequency";
          "samplinginstants";
          "domain";
          "timedomain"};

    keys = [keys; ka];
  endif

endfunction
