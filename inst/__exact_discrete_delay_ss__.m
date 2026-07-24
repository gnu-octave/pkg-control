## -*- texinfo -*-
## Internal helper: exact causal discrete-time ss realization of z^-k (a
## pure k-sample delay chain). Used by pade() (discrete-time InternalDelay
## closure) and c2d's DelayModeling="state" InternalDelay handling -- both
## need an EXACT (not approximate) rational discrete delay, since integer
## sample delays are already exactly rational; unlike the continuous case,
## no Pade approximation is needed or appropriate here.
function dsys = __exact_discrete_delay_ss__ (k, tsam)

  if (k == 0)
    dsys = ss (zeros (0, 0), zeros (0, 1), zeros (1, 0), 1, tsam);
    return;
  endif

  A = diag (ones (k-1, 1), 1);
  B = [zeros(k-1, 1); 1];
  C = [1, zeros(1, k-1)];
  D = 0;

  dsys = ss (A, B, C, D, tsam);

endfunction

%!test  # k=0: pure passthrough, no states
%! d = __exact_discrete_delay_ss__ (0, 0.1);
%! assert (isa (d, "ss"));
%! assert (isdt (d));
%! [a, b, c, dd] = ssdata (d);
%! assert (rows (a), 0);
%! assert (dd, 1);

%!test  # k=3: exact z^-3, verified via freqresp against the analytic z^-3 phase/magnitude
%! tsam = 0.1;
%! d = __exact_discrete_delay_ss__ (3, tsam);
%! assert (isa (d, "ss"));
%! [a, b, c, dd] = ssdata (d);
%! assert (rows (a), 3);
%! w = [0.1, 1, 5];
%! z = exp (1i * w * tsam);
%! expected = z .^ (-3);
%! assert (reshape (freqresp (d, w), 1, []), expected, 1e-10);

%!test  # step response: output should be exactly the (delayed-by-3) unit step
%! tsam = 0.1;
%! d = __exact_discrete_delay_ss__ (3, tsam);
%! N = 10;
%! u = ones (N, 1);
%! t = (0:N-1)' * tsam;
%! y = lsim (d, u, t);
%! expected = [zeros(3,1); ones(N-3,1)];
%! assert (y, expected, 1e-10);
