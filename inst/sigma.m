## Copyright (C) 2009 Lukas Reichlin <lukas.reichlin@swissonline.ch>
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn{Function File} {[@var{sv}, @var{w}] =} sigma (@var{sys})
## @deftypefnx{Function File} {[@var{sv}, @var{w}] =} sigma (@var{sys}, @var{w})
## @deftypefnx{Function File} {[@var{sv}, @var{w}] =} sigma (@var{sys}, @var{[]}, @var{ptype})
## @deftypefnx{Function File} {[@var{sv}, @var{w}] =} sigma (@var{sys}, @var{w}, @var{ptype})
## If no output arguments are given, the singular value plot of a MIMO
## system is printed on the screen;
## otherwise, the singular values of the system data structure are
## computed and returned.
##
## @strong{Inputs}
## @table @var
## @item sys
## System data structure. Must be either purely continuous or discrete;
## see @code{is_digital}.
## @item w
## Optional vector of frequency values. If @var{w} is not specified, it
## will be calculated by @code{bode_bounds}.
## @item ptype = 0
## Singular values of the frequency response H of system sys. Default Value.
## @item ptype = 1
## Singular values of the frequency response inv(H); i.e. inversed system.
## @item ptype = 2
## Singular values of the frequency response I + H; i.e. inversed sensitivity
## (or return difference) if H = P * C.
## @item ptype = 3
## Singular values of the frequency response I + inv(H); i.e. inversed complementary
## sensitivity if H = P * C.
## @end table
##
## @strong{Outputs}
## @table @var
## @item sv
## Array of singular values. For a system with m inputs and p outputs, the array sv
## has min(m,p) rows and as many columns as frequency points (length of w).
## The singular values at the frequency w(k) are given by sv(:,k).
## @item w
## Vector of frequency values used.
## @end table
##
## @seealso{bode, svd, bode_bounds, is_digital}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@swissonline.ch>
## Version: 0.1.2

function [sv_r, w_r] = sigma (sys, w = [], ptype = 0)

  ## Check whether arguments are OK
  if (nargin < 1 || nargin > 3)
    print_usage ();
  endif

  if(! isstruct (sys))
    error ("sigma: first argument sys must be a system data structure");
  endif

  if (! isvector (w) && ! isempty (w))
    error ("sigma: second argument w must be a vector of frequencies");
  endif

  if (is_siso (sys))
    warning ("sigma: sys is a SISO system. You might want to use bode(sys) instead.");
  endif

  ## Get state space matrices
  sys = sysupdate (sys, "ss");
  [A, B, C, D] = sys2ss (sys);
  I = eye (size (A));

  ## Get system information
  [n_c, n_d, m, p] = sysdimensions (sys);
  digital = is_digital (sys, 2);
  t_sam = sysgettsam (sys);

  ## Error for mixed systems
  if (digital == -1)
    error ("sigma: system must be either purely continuous or purely discrete");
  endif

  ## Handle plot type
  if (nargin == 3)
    if (isfloat (ptype))  # Numeric constants like 2 are NOT integers in Octave!
      if (ptype != 0 && ptype != 1 && ptype != 2 && ptype != 3)
        error ("sigma: third argument ptype must be 0, 1, 2 or 3");
      endif
    else
      error ("sigma: third argument ptype must be a number");
    endif

    if (m != p)
      error ("sigma: system must be square for ptype 1, 2 or 3");
    endif
    J = eye (m);
  endif

  ## Find interesting frequency range w if not specified
  if (isempty (w))

    ## Since no frequency vector w has been specified, the interesting
    ## decades of the sigma plot have to be found. The already existing
    ## function bode_bounds does exactly that.

    zer = tzero (sys);
    pol = eig (A);

    ## Begin plot at 10^dec_min, end plot at 10^dec_max [rad/s]
    [dec_min, dec_max] = bode_bounds (zer, pol, digital, t_sam);

    w = logspace (dec_min, dec_max, 1000);  # [rad/s]
  endif

  if (digital)  # Discrete system
    s = exp (i * w * t_sam);
  else  # Continuous system
    s = i * w;
  endif

  l_s = length (s);
  sv = zeros (min (m, p), l_s);

  switch (ptype)
    case 0  # Default system
      for k = 1 : l_s  # Repeat for every frequency s
        H = C * inv (s(k)*I - A) * B  +  D;  # Frequency Response Matrix
        sv(:, k) = svd (H);  # Singular Value Decomposition
      endfor

    case 1  # Inversed system
      for k = 1 : l_s
        H = inv (C * inv (s(k)*I - A) * B  +  D);
        sv(:, k) = svd (H);
      endfor

    case 2  # Inversed sensitivity
      for k = 1 : l_s
        H = J  +  C * inv (s(k)*I - A) * B  +  D;
        sv(:, k) = svd (H);
      endfor

    case 3  # Inversed complementary sensitivity
      for k = 1 : l_s
        H = J  +  inv (C * inv (s(k)*I - A) * B  +  D);
        sv(:, k) = svd (H);
      endfor
  endswitch

  if (nargout == 0)  # Plot the information

    ## Convert to dB for plotting
    sv_db = 20 * log10 (sv);

    ## Determine axes
    ax_vec = __axis2dlim__ ([[w(:), min(sv_db)(:)]; [w(:), max(sv_db)(:)]]);
    ax_vec(1:2) = [min(w), max(w)];

    ## Determine xlabel
    if (digital)
      xl_str = sprintf ("Frequency [rad/s]     Pi / T = %g", pi/t_sam);
    else
      xl_str = "Frequency [rad/s]";
    endif

    ## Plot results
    semilogx (w, sv_db, "b")
    axis (ax_vec)
    title ("Singular Values")
    xlabel (xl_str)
    ylabel ("Singular Values [dB]")
    grid on
  else  # Return values
    sv_r = sv;
    w_r = w;
  endif

endfunction


%!shared sv_exp, w_exp, sv_obs, w_obs
%! A = [1 2; 3 4];
%! B = [5 6; 7 8];
%! C = [4 3; 2 1];
%! D = [8 7; 6 5];
%! w = [2 3];
%! sv_exp = [7.91760889977901, 8.62745836756994; 0.698526948925716, 0.608629874340667];
%! w_exp = [2 3];
%! [sv_obs, w_obs] = sigma (ss (A, B, C, D), w);
%!assert (sv_obs, sv_exp, 16*eps); # tolerance manually tweaked
%!assert (w_obs, w_exp, 2*eps);
