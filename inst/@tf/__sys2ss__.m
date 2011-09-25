## Copyright (C) 2009, 2011   Lukas F. Reichlin
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
## TF to SS conversion.

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.2

function [retsys, retlti] = __sys2ss__ (sys)

  ## TODO: determine appropriate tolerance from number of inputs
  ##       (since we multiply all denominators in a row), index, ...
  ##       default tolerance from TB01UD is TOLDEF = N*N*EPS 

  ## SECRET WISH: a routine which accepts individual denominators for
  ##              each channel and which supports descriptor systems

  [p, m] = size (sys);
  [num, den] = tfdata (sys);

  ## new cells for the TF of same row denominators
  numc = cell (p, m);
  denc = cell (p, 1);
  
  ## multiply all denominators in a row and
  ## update each numerator accordingly 
  for i = 1 : p
    denc(i) = __conv__ (den{i,:});
    for j = 1 : m
      idx = setdiff (1:m, j);
      numc(i,j) = __conv__ (num{i,j}, den{i,idx});
    endfor
  endfor

  len_numc = cellfun (@length, numc);
  len_denc = cellfun (@length, denc);

  ## check for properness  
  ## tfpoly ensures that there are no leading zeros
  tmp = len_numc > repmat (len_denc, 1, m);
  if (any (tmp(:)))
    error ("tf: tf2ss: system must be proper");
  endif

  ## create arrays and fill in the data
  ## in a way that Slicot TD04AD can use
  max_len_denc = max (len_denc(:));
  ucoeff = zeros (p, m, max_len_denc);
  dcoeff = zeros (p, max_len_denc);
  index = len_denc-1;

  for i = 1 : p
    len = len_denc(i);
    dcoeff(i, 1:len) = denc{i};
    for j = 1 : m
      ucoeff(i, j, len-len_numc(i,j)+1 : len) = numc{i,j};
    endfor
  endfor

  [a, b, c, d] = sltd04ad (ucoeff, dcoeff, index, sqrt (eps));

  retsys = ss (a, b, c, d);
  retlti = sys.lti;   # preserve lti properties such as tsam

endfunction


## convolution for more than two arguments
function vec = __conv__ (vec, varargin)

  if (nargin == 1)
    return;
  else
    for k = 1 : nargin-1
      vec = conv (vec, varargin{k});
    endfor
  endif

endfunction
