## Copyright (C) 2013   Lukas F. Reichlin
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
## @deftypefn{Function File} {@var{P} =} mktito (@var{P}, @var{nmeas}, @var{ncon})
## Partition @acronym{LTI} plant @var{P} for robust controller synthesis.
##
## @strong{Inputs}
## @table @var
## @item P
## Generalized plant.
## @item nmeas
## Number of measured outputs v.  The last @var{nmeas} outputs of @var{P} are connected to the
## inputs of controller @var{K}.  The remaining outputs z (indices 1 to p-nmeas) are used
## to calculate the H-2/H-infinity norm.
## @item ncon
## Number of controlled inputs u.  The last @var{ncon} inputs of @var{P} are connected to the
## outputs of controller @var{K}.  The remaining inputs w (indices 1 to m-ncon) are excited
## by a harmonic test signal.
## @end table
##
## @strong{Outputs}
## @table @var
## @item P
## Partitioned plant.
## @end table
##
## @strong{Block Diagram}
## @example
## @group
##
## min||N(K)||             N = lft (P, K)
##  K         norm
##
##                +--------+  
##        w ----->|        |-----> z
##                |  P(s)  |
##        u +---->|        |-----+ v
##          |     +--------+     |
##          |                    |
##          |     +--------+     |
##          +-----|  K(s)  |<----+
##                +--------+
##
##                +--------+      
##        w ----->|  N(s)  |-----> z
##                +--------+
## @end group
## @end example
##
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2013
## Version: 0.1

function P = mktito (P, nmeas, ncon)

  if (nargin != 3)
    print_usage ();
  endif

  if (! isa (P, "lti"))
    error ("mktito: first argument must be an LTI model");
  endif
  
  [p, m] = size (P);
  
  ## TODO: improve argument checking
  if (! is_real_scalar (nmeas))
    error ("mktito: second argument invalid");
  endif
  
  if (! is_real_scalar (ncon))
    error ("mktito: third argument invalid");
  endif
  
  outgroup = struct ("Z", 1:p-nmeas, "V", p-nmeas+1:p);
  outname = vertcat (strseq ("z", 1:p-nmeas), strseq ("v", 1:nmeas));
  
  ingroup = struct ("W", 1:m-ncon, "U", m-ncon+1:m);
  inname = vertcat (strseq ("w", 1:m-ncon), strseq ("u", 1:ncon));

  P = set (P, "outgroup", outgroup, "ingroup", ingroup, ...
              "outname", outname, "inname", inname);

endfunction
