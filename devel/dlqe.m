function [l, p, z, e] = dlqe (a, g, c, q, r, s = [])

  if (nargin < 5 || nargin > 6)
    print_usage ();
  endif

  if (isempty (s))
    [~, p, e] = dlqr (a.', c.', g*q*g.', r);
  else
    [~, p, e] = dlqr (a.', c.', g*q*g.', r, g*s);
  endif

  ## k computed by dlqr would be
  ## k = (r + c*p*c.') \ (c*p*a.' + s.')
  ## such that  l = a \ k.'
  ## what about the s term?
    
  l = p*c.' / (c*p*c.' + r);
  ## z = p - p*c.' / (c*p*c.' + r) * c*p;
  z = p - l*c*p;

endfunction