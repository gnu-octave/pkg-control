function [retsys, lscale, rscale] = __prescale__ (sys, optarg = 0.0)

  if (isempty (sys.e))
    [a, b, c, ~, scale] = sltb01id (sys.a, sys.b, sys.c, optarg);
    retsys = ss (a, b, c, sys.d);
    lscale = rscale = scale;
  else
    [a, e, b, c, lscale, rscale] = sltg01ad (sys.a, sys.e, sys.b, sys.c, optarg);
    retsys = dss (a, b, c, sys.d, e);
  endif
  
  retsys.scaled = true;
  retsys.lti = sys.lti;  # retain i/o names and tsam

endfunction
