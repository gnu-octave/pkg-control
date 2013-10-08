## Copyright (C) 2009, 2013   Lukas F. Reichlin
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
## @deftypefn {Function File} {@var{sys} =} connect (@var{sys}, @var{cm}, @var{inputs}, @var{outputs})
## @deftypefnx {Function File} {@var{sys} =} connect (@var{sys1}, @var{sys2}, @dots{}, @var{sysN}, @var{inputs}, @var{outputs})
## Arbitrary interconnections between the inputs and outputs of an @acronym{LTI} model.
## @seealso{sumblk}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2009
## Version: 0.2

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

    ## TODO: proper argument checking
    ## TODO: replace nested for-if statement

    ## if (! is_real_matrix (cm))
    ##   error ("connect: second argument must be a matrix with real-valued coefficients");
    ## endif

    M = zeros (m, p);
    in = cm(:, 1);
    out = cm(:, 2:cmcols);

    for a = 1 : cmrows
      for b = 1 : cmcols-1
        if (out(a, b) != 0)
          M(in(a, 1), abs (out(a, b))) = sign (out(a, b));
        endif
      endfor
    endfor

  else                                  # connect (sys1, sys2, ..., sysN, in_idx, out_idx)

    lti_idx = cellfun (@isa, varargin, {"lti"});
    sys = append (varargin{lti_idx});
    io_idx = ! lti_idx;
    
    if (nnz (io_idx) != 2)
      error ("connect: require arguments 'inputs' and 'outputs'");
    endif
    
    in_idx = varargin(io_idx){1};
    out_idx = varargin(io_idx){2};

    inname = sys.inname;
    outname = sys.outname;
    
    ioname = intersect (inname, outname);
    
    tmp = cellfun (@(x) find (strcmp (inname, x)(:)), ioname, "uniformoutput", false);
    inputs = vertcat (tmp{:});  # there could be more than one input with the same name
    
    ## FIXME: sys_prune will error out if names in out_idx and in_idx are not unique
    ##        the dark side handles cases with common in_idx names
    
    [p, m] = size (sys);
    M = zeros (m, p);
    for k = 1 : length (inputs)
      outputs = strcmp (outname, inname(inputs(k)));
      M(inputs(k), :) = outputs;
    endfor

  endif

  sys = __sys_connect__ (sys, M);
  sys = __sys_prune__ (sys, out_idx, in_idx);

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
