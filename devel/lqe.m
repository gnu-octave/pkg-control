function [g, x, l] = lqe (a, c, q, r = [], s = [], e = [])

  if (nargin < 3 || nargin > 6)
    print_usage ();
  endif

  if (isa (a, "lti"))
    [g, x, l] = lqr (a.', c, q, r);  # lqr (sys.', q, r, s)
  else
    [g, x, l] = lqr (a.', c.', q, r, s, e.');
  endif
  
  g = g.'

endfunction