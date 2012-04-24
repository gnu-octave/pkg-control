function [l, p, e] = lqe (a, g, c, q = [], r = [], s = [])

  if (nargin < 3 || nargin > 6)
    print_usage ();
  endif

  if (isa (a, "lti"))
    [l, p, e] = lqr (a.', g, c, q);         # lqe (sys.', q, r, s), g=I, works like  lqr (sys.', q, r, s).'
  elseif (isempty (g))
    [l, p, e] = lqr (a.', c.', q, r, s);    # lqe (a, [], c, q, r, s), g=I, works like  lqr (a.', c.', q, r, s).'
  elseif (isempty (s))
    [l, p, e] = lqr (a.', c.', g*q*g.', r);
  else
    [l, p, e] = lqr (a.', c.', g*q*g.', r, g*s);
  endif

  l = l.';

endfunction