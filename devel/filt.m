function sys = filt (num = {}, den = {}, tsam = -1, varargin)

  if (! iscell (num))
    num = {num};
  endif
  
  if (! iscell (den))
    den = {den};
  endif

  ## shall i remove trailing zeros before postpadding?
  num = cellfun (@__remove_trailing_zeros__, num, "uniformoutput", false);
  den = cellfun (@__remove_trailing_zeros__, den, "uniformoutput", false);

  lnum = cellfun (@length, num, "uniformoutput", false);
  lden = cellfun (@length, den, "uniformoutput", false);

  lmax = cellfun (@max, lnum, lden, "uniformoutput", false);

  num = cellfun (@postpad, num, lmax, "uniformoutput", false);
  den = cellfun (@postpad, den, lmax, "uniformoutput", false);
  
  sys = tf (num, den, tsam, varargin{:});

endfunction


function p = __remove_trailing_zeros__ (p)

  idx = find (p != 0);
  
  if (isempty (idx))
    p = 0;
  else
    p = p(1 : idx(end));
  endif

endfunction