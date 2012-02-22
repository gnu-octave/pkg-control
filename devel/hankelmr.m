## Copyright (C) 2012   Lukas F. Reichlin
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
## @deftypefn{Function File} {[@var{Gr}, @var{info}] =} hankelmr (@var{G}, @dots{})
## Wrapper for @command{hnamodred}.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: January 2012
## Version: 0.1

function [Gr, info] = hankelmr (varargin)

  persistent warned = false;
  if (! warned)
    warned = true;
    warning ("control:wrapper",
             "control: 'hankelmr' is just a compatibility wrapper for 'hnamodred'");
  endif

  [Gr, info] = hnamodred (varargin{:});

endfunction
