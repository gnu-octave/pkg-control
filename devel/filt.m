function sys = filt (num = {}, den = {}, tsam = -1)

  lnum = cellfun (@length, num, "uniformoutput", false);
  lden = cellfun (@length, den, "uniformoutput", false);

  lmax = cellfun (@max, lnum, lden, "uniformoutput", false);

  num = cellfun (@postpad, num, lmax, "uniformoutput", false);
  den = cellfun (@postpad, den, lmax, "uniformoutput", false);
  
  sys = tf (num, den, tsam);

endfunction