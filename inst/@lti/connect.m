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
## @deftypefn {Function File} {@var{sys} =} connect (@var{sys1}, @var{sys2}, @dots{}, @var{sysN}, @var{inputs}, @var{outputs})
## @deftypefnx {Function File} {@var{sys} =} connect (@var{sys}, @var{cm}, @var{inputs}, @var{outputs})
## Name-based or index-based interconnections between the inputs and outputs of @acronym{LTI} models.
## 
## @strong{Inputs}
## @table @var
## @item sys1, @dots{}, sysN
## @acronym{LTI} models to be connected.  The properties 'inname' and 'outname'
## of each model should be set according to the desired input-output connections.
## @item inputs
## For name-based interconnections, string or cell of strings containing the names
## of the inputs to be kept.  The names must be part of the properties 'ingroup' or
## 'inname'.  For index-based interconnections, vector containing the indices of the
## inputs to be kept.
## @item outputs
## For name-based interconnections, string or cell of strings containing the names
## of the outputs to be kept.  The names must be part of the properties 'outgroup' 
## or 'outname'.  For index-based interconnections, vector containing the indices of
## the outputs to be kept.
## @item cm
## Connection matrix (not name-based).  Each row of the matrix represents a summing
## junction.  The first column holds the indices of the inputs to be summed with
## outputs of the subsequent columns.  The output indices can be negative, if the output
## is to be substracted, or zero.  For example, the row
## @example
## [2 0 3 -4 0]
## @end example
## or
## @example
## [2 -4 3]
## @end example
## will sum input u(2) with outputs y(3) and y(4) as
## @example
## u(2) + y(3) - y(4).
## @end example
## @end table
##
## @strong{Outputs}
## @table @var
## @item sys
## Resulting interconnected system with outputs @var{outputs} and
## inputs @var{inputs}.
## @end table
##
## @seealso{sumblk}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.4

function sys = connect (varargin)

  if (nargin < 2)
    print_usage ();
  endif
  
  if (is_real_matrix (varargin{2}))     # connect (sys, cm, in_idx, out_idx)
  
    if (nargin != 4)
      print_usage ();
    endif
    
    sys = varargin{1};
    cm = varargin{2};
    in_idx = varargin{3};
    out_idx = varargin{4};
    
    [p, m] = size (sys);
    [cmrows, cmcols] = size (cm);

    ## if (! is_real_matrix (cm))
    ##   error ("connect: second argument must be a matrix with real-valued coefficients");
    ## endif

    M = zeros (m, p);
    in = cm(:, 1);
    out = cm(:, 2:cmcols);
    
    ## check sizes and integer values
    if (! isequal (cm, floor (cm)))
      error ("connect: matrix 'cm' must contain integer values (index-based interconnection)");
    endif
    
    if ((min (in) <= 0) || (max (in) > m))
      error ("connect: 'cm' input index in out of range (index-based interconnection)");
    endif
    
    if (max (abs (out(:))) > p)
      error ("connect: 'cm' output index out of range (index-based interconnection)");
    endif
    
    if ((! is_real_vector (in_idx)) || (! isequal (in_idx, floor (in_idx))))
      error ("connect: 'inputs' must be a vector of integer values (index-based interconnection)");
    endif
    
    if ((max (in_idx) > m) || (min (in_idx) <= 0))
      error ("connect: index in vector 'inputs' out of range (index-based interconnection)");
    endif
    
    if ((! is_real_vector (out_idx)) || (! isequal (out_idx, floor (out_idx))))
      error ("connect: 'outputs' must be a vector of integer values (index-based interconnection)");
    endif
    
    if ((max (out_idx) > p) || (min (out_idx) <= 0))
      error ("connect: index in vector 'outputs' out of range (index-based interconnection)");
    endif
    
    for a = 1 : cmrows
      out_tmp = out(a, (out(a,:) != 0));
      if (! isempty (out_tmp))
        M(in(a,1), abs (out_tmp)) = sign (out_tmp);
      endif
    endfor

    sys = __sys_connect__ (sys, M);
    sys = __sys_prune__ (sys, out_idx, in_idx);

  else                                  # connect (sys1, sys2, ..., sysN, in_idx, out_idx)

    lti_idx = cellfun (@isa, varargin, {"lti"});
    sys = blkdiag (varargin{lti_idx});
    io_idx = ! lti_idx;
    
    if (nnz (io_idx) == 2)
      in_idx = varargin(io_idx){1};
      out_idx = varargin(io_idx){2};      
    else
      in_idx = ":";
      out_idx = ":";
    endif

    inname = sys.inname;
    if (any (cellfun (@isempty, inname)))
      error ("connect: all inputs must have names");
    endif

    outname = sys.outname;
    if (any (cellfun (@isempty , outname)))
      error ("connect: all outputs must have names");
    endif
    
    ioname = intersect (inname, outname);
    
    tmp = cellfun (@(x) find (strcmp (inname, x)(:)), ioname, "uniformoutput", false);
    inputs = vertcat (tmp{:});  # there could be more than one input with the same name
    
    [p, m] = size (sys);
    M = zeros (m, p);
    for k = 1 : length (inputs)
      outputs = strcmp (outname, inname(inputs(k)));
      M(inputs(k), :) = outputs;
    endfor

    sys = __sys_connect__ (sys, M);

    ## sys_prune will error out if names in out_idx and in_idx are not unique
    ## the dark side handles cases with common in_idx names - so do we

    inname_u = unique (inname);
    if (numel (inname_u) != numel (inname))
      tmp = cellfun (@(u) strcmp (u, inname), inname_u, "uniformoutput", false);
      mat = double (horzcat (tmp{:}));
      scl = ss (mat, "inname", inname_u, "outname", inname);
      sys = sys * scl;
      if (is_real_vector (in_idx))
        warning ("connect: use names instead of indices for argument 'inputs'");
      endif
    endif

    sys = __sys_prune__ (sys, out_idx, in_idx);

    if (isa (sys, "ss"))
      sys = sminreal (sys);
    endif

  endif

endfunction


%!shared T, Texp
%! P = Boeing707;
%! I = ss (-eye (2));
%! I.inname = P.outname;
%! I.outname = P.inname;
%! T = connect (P, I, P.inname, P.outname);
%! Texp = feedback (P);
%!assert (T.a, Texp.a, 1e-4);
%!assert (T.b, Texp.b, 1e-4);
%!assert (T.c, Texp.c, 1e-4);
%!assert (T.d, Texp.d, 1e-4);
