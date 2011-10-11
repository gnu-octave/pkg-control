## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: October 2011
## Version: 0.1

function [p, m] = __iddata_dim__ (y, u)

  if (! is_real_matrix (y, u))
    error ("iddata: inputs and outputs must be real");
  endif
  
  [ly, p] = size (y);
  [lu, m] = size (u);
  
  if (ly != lu)
    error ("iddata: matrices ""y"" and ""u"" must have the same number of samples (rows)");
  endif

  if (ly < p)
    warning ("iddata: more outputs than samples - matrice ""y"" should probably be transposed");
  endif
  
  if (lu < m)
    warning ("iddata: more inputs than samples - matrice ""u"" should probably be transposed");
  endif

endfunction