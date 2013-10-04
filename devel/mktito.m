function sys = mktito (sys, nmeas, ncont)

  if (nargin != 3)
    print_usage ();
  endif

  if (! isa (sys, "lti"))
    error ("mktito: first argument must be an LTI model");
  endif
  
  [p, m] = size (sys);
  
  ## TODO: improve argument checking
  if (! is_real_scalar (nmeas, ncont))
    print_usage ();
  endif
  
  outgroup = struct ("Y1", 1:p-nmeas, "Y2", p-nmeas+1:p);
  ingroup = struct ("U1", 1:m-ncont, "U2", m-ncont+1:m);

  sys = set (sys, "outgroup", outgroup, "ingroup", ingroup);

endfunction
