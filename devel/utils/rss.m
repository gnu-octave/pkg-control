function sys = rss (n, p = 1, m = 1)

  if (nargin < 1 || nargin > 3 || ! is_real_scalar (n, p, m))
    print_usage ();
  endif

  a = rand (n, n);
  a = 0.5 * (a + a.');
  a(1:n+1:end) = -n;
  b = randn (n, m);
  c = randn (p, n);
  d = randn (p, m);

  sys = ss (a, b, c, d);

endfunction
