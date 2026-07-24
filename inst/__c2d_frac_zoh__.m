## -*- texinfo -*-
## Internal helper: exact zero-order-hold discretization of one input
## channel's fractional (non-integer-sample) continuous delay.
##
## Given continuous state matrix A, one column b_col of the continuous B
## matrix, sample time tsam, and that channel's total continuous delay tau
## (seconds), returns the discrete Ad = expm(A*tsam), the discrete input
## column Bd to use together with the ordinary d = floor(tau/tsam) sample
## delay (already handled by the existing InputDelay/IODelay machinery), and
## an `extra` struct describing the one extra register state needed to
## carry the sub-sample remainder exactly (Astrom & Wittenmark's split-ZOH
## construction), or extra.needed = false when tau is already an integer
## number of samples.
function [Ad, Bd, extra] = __c2d_frac_zoh__ (A, b_col, tsam, tau)

  d = floor (tau / tsam + 1e-9);      # guard against tau/tsam landing just
                                       # under an integer due to roundoff
  tau_frac = tau - d * tsam;

  [Ad, G_full] = __sl_mb05nd__ (A, tsam, eps);

  if (tau_frac < 1e-8 * tsam)
    Bd = G_full * b_col;
    extra = struct ("needed", false, "d", d, "g0", []);
    return;
  endif

  [Ad_short, G1] = __sl_mb05nd__ (A, tsam - tau_frac, eps);
  [~, Gf] = __sl_mb05nd__ (A, tau_frac, eps);

  Bd = G1 * b_col;
  g0 = Ad_short * Gf * b_col;

  extra = struct ("needed", true, "d", d, "g0", g0);

endfunction

%!test  # scalar plant: verify against direct numerical integration
%! A = -2; b = 1; tsam = 0.1; tau = 0.25;
%! [Ad, Bd, extra] = __c2d_frac_zoh__ (A, b, tsam, tau);
%! assert (extra.needed, true);
%! assert (extra.d, 2);
%! # brute-force check: Ad == expm(A*tsam)
%! assert (Ad, expm (A*tsam), 1e-10);
%! # brute-force check of Bd (=G1*b) and extra.g0 (=G0*b) via direct
%! # numerical quadrature of the defining integrals
%! tau_frac = tau - extra.d*tsam;
%! G1_direct = integral (@(s) expm (A*(tsam-tau_frac-s)), 0, tsam-tau_frac, "ArrayValued", true);
%! G0_direct = integral (@(s) expm (A*s), tsam-tau_frac, tsam, "ArrayValued", true);
%! assert (Bd, G1_direct*b, 1e-8);
%! assert (extra.g0, G0_direct*b, 1e-8);
%!
%!test  # zero fractional part: extra.needed is false, Bd matches plain zoh
%! A = -2; b = 1; tsam = 0.1; tau = 0.3;      # exactly 3 samples
%! [Ad, Bd, extra] = __c2d_frac_zoh__ (A, b, tsam, tau);
%! assert (extra.needed, false);
%! assert (extra.d, 3);
%! [~, G] = __sl_mb05nd__ (A, tsam, eps);
%! assert (Bd, G*b, 1e-10);
