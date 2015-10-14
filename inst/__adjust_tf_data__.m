function [num, den, tsam, tfvar] = __adjust_tf_data__ (num, den, tsam)

  if (is_real_matrix (num) && isempty (den))    # static gain  tf (4),  tf (matrix)
    num = num2cell (num);
    tsam = -2;
  endif

  if (! iscell (num))
    num = {num};
  endif
  
  if (! iscell (den))
    den = {den};
  endif

  num = cellfun (@tfpoly, num, "uniformoutput", false);
  den = cellfun (@tfpoly, den, "uniformoutput", false);
  
  if (tsam == 0)
    tfvar = "s";
  elseif (tsam == -2)
    tfvar = "x";
  else
    tfvar = "z";
  endif

endfunction
