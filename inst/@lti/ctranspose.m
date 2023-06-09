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
## @deftypefn{Function File} {@var{SYST} =} ctranspose (@var{SYS})
## Conjugate transpose or pertransposition of @acronym{LTI} objects.
## Used by Octave for "sys'".
## For a transfer-function matrix G, G' denotes the conjugate
## of G given by G.'(-s) for a continuous-time system or G.'(1/z)
## for a discrete-time system.
## The frequency response of the pertransposition of G is the
## Hermitian (conjugate) transpose of G(jw), i.e.
## freqresp (G', w) = freqresp (G, w)'.
## @strong{WARNING:} Do @strong{NOT} use this for dual problems,
## use the transpose "sys.'" (note the dot) instead.
##
## @strong{Inputs}
## @table @var
## @item SYS
## System to be transposed.
## @end table
##
## @strong{Outputs}
## @table @var
## @item SYST
## Conjugate transposed of @var{SYS}.
## @end table
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: May 2012
## Version: 0.2

function sys = ctranspose (sys)

  if (nargin != 1)  # prevent sys = ctranspose (sys1, sys2, sys3, ...)
    error ("lti: ctranspose: this is an unary operator");
  endif

  [p, m] = size (sys);
  ct = isct (sys);

  sys = __ctranspose__ (sys, ct);

  sys.inname = repmat ({""}, p, 1);
  sys.outname = repmat ({""}, m, 1);
  sys.ingroup = struct ();
  sys.outgroup = struct ();

endfunction
