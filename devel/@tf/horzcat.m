function sys = horzcat (sys, varargin)

  warning ("tf: horzcat!");
  
  sys = tf (sys);
  varargin = cellfun (@tf, varargin, "uniformoutput", false);

  for k = 1 : (nargin-1)
  
    sys1 = sys;
    sys2 = varargin{k};
    
    sys = tf ();
    sys.lti = __lti_group__ (sys1.lti, sys2.lti, "horz");
    
    [p1, m1] = size (sys1.num);
    [p2, m2] = size (sys2.num);
    
    if (p1 != p2)
      error ("tf: horzcat: number of system outputs incompatible: [(%dx%d), (%dx%d)]",
              p1, m1, p2, m2);
    endif
    
    sys.num = [sys1.num, sys2.num];
    sys.den = [sys1.den, sys2.den];
    
    if (sys1.tfvar == sys2.tfvar)
      sys.tfvar = sys1.tfvar;
    elseif (sys1.tfvar == "x")
      sys.tfvar = sys2.tfvar;
    else
      sys.tfvar = sys1.tfvar;
    endif

    if (sys1.inv || sys2.inv)
      sys.inv = true;
    endif

  endfor

endfunction
