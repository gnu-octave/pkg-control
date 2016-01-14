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

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2011
## Version: 0.1

function [p, m, e] = __iddata_dim__ (y, u)

  e = numel (y);                        # number of experiments

  if (isempty (u))                      # time series data, outputs only
    [p, m] = cellfun (@__experiment_dim__, y, "uniformoutput", false);
  elseif (e == numel (u))               # outputs and inputs present
    [p, m] = cellfun (@__experiment_dim__, y, u, "uniformoutput", false);
  else
    error ("iddata: require input and output data with matching number of experiments");
  endif

  if (e > 1 && ! isequal (p{:}))
    error ("iddata: require identical number of output channels for all experiments");
  endif

  if (e > 1 && ! isequal (m{:}))
    error ("iddata: require identical number of input channels for all experiments");
  endif

  p = p{1};
  m = m{1};

endfunction


function [p, m] = __experiment_dim__ (y, u = [])

  if (! is_matrix (y, u))
    error ("iddata: inputs and outputs must be real or complex matrices");
  endif
  
  [ly, p] = size (y);
  [lu, m] = size (u);
  
  if (! isempty (u) && ly != lu)
    error ("iddata: matrices 'y' (%dx%d) and 'u' (%dx%d) must have the same number of samples (rows)", ...
           ly, p, lu, m);
  endif

  if (ly < p)
    warning ("iddata:transpose", "iddata: more outputs than samples - matrice 'y' should probably be transposed");
  endif
  
  if (lu < m)
    warning ("iddata:transpose", "iddata: more inputs than samples - matrice 'u' should probably be transposed");
  endif

endfunction
