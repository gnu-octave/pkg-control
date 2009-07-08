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
## @deftypefn{Function File} {[@var{sigma_min}, @var{sigma_max}, @var{w}] =} svplot (@var{sys})
## @deftypefnx{Function File} {[@var{sigma_min}, @var{sigma_max}, @var{w}] =} svplot (@var{sys}, @var{w})
## If no output arguments are given, the singular value plot of a MIMO
## system over a range of frequencies is printed on the screen;
## otherwise, the singular values of the system data structure are
## computed and returned.
##
## @strong{Inputs}
## @table @var
## @item sys
## System data structure (must be either purely continuous or discrete;
## see @code{is_digital}).
## @item w
## Optional vector of frequency values. If @var{w} is not specified, it
## will be calculated by @code{bode_bounds}.
## @end table
##
## @strong{Outputs}
## @table @var
## @item sigma_min
## Vector of minimal singular values.
## @item sigma_max
## Vector of maximal singular values.
## @item w
## Vector of frequency values used.
## @end table
##
## @seealso{is_digital}
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@swissonline.ch>
## Version: 0.1


function [sigma_min_r, sigma_max_r, w_r] = svplot (sys, w)

  ## Check whether arguments are OK
  if (nargin < 1 || nargin > 2)
    print_usage ();
  endif

  if(! isstruct (sys))
    error ("first argument sys must be a system data structure");
  endif

  if (nargin == 2)
    if (! isvector (w))
      error ("second argument w must be a vector of frequencies");
    endif
  endif

  ## Get State Space Matrices
  sys = sysupdate (sys, "ss");
  [A, B, C, D] = sys2ss (sys);
  I = eye (size (A));
  
  ## Get system information
  [n_c_states, n_d_states, n_in, n_out] = sysdimensions (sys);
  
  ## Warning for discrete systems, see FIXME below!
  if (n_d_states != 0)
    warning ("discrete systems have not been tested!");
  endif

  ## Find interesting frequency range w if not specified
  if (nargin == 1)

    ## Since no frequency vector w has been specified, the interesting
    ## decades of the sigma plot have to be found. The already existing
    ## function bode_bounds does exactly that, unfortunately for SISO
    ## systems only. Therefore, a SISO system from every input m to
    ## every output p is created. Then for every SISO system the
    ## interesting frequency range is calculated by bode_bounds.
    ## Afterwards, the min/max decade can be found by the min()/max()
    ## command.

    j = 1; # Index counter (number of SISO systems)
    
    for m = 1 : n_in # For all inputs
      for p = 1 : n_out # For all outputs

        sysp = sysprune (sys, p, m); # Create SISO system
        [zer, pol, k, tsam, inname, outname] = sys2zp (sysp);

        if (tsam == 0) # Continuous system
          DIGITAL = 0;
        else # Discrete system
          DIGITAL = 1;
          ## FIXME: Discrete systems not yet tested!
        endif

        [wmin, wmax] = bode_bounds (zer, pol, DIGITAL, tsam);

        ## Save result for later evaluation
        w_min(j) = wmin; 
        w_max(j) = wmax;
        j++;
      endfor
    endfor

    dec_min = min (w_min); # Begin plot at 10^dec_min [rad/s]
    dec_max = max (w_max); # End plot at 10^dec_min [rad/s]
    n_steps = 1000; # Number of frequencies evaluated for plotting

    w = logspace (dec_min, dec_max, n_steps); # [rad/s]
  endif


  ## Repeat for every frequency w
  for k = 1 : size (w, 2)

    ## Frequency Response Matrix
    P = C * inv (i*w(k)*I - A) * B  +  D;

    ## Singular Value Decomposition
    sigma = svd (P);
    sigma_min(k) = min (sigma);
    sigma_max(k) = max (sigma);
  endfor

  if (nargout == 0) # Plot the information

    ## Convert to dB for plotting
    sigma_min_db = 20 * log10 (sigma_min);
    sigma_max_db = 20 * log10 (sigma_max);
   
    ## Plot results
    semilogx (w, sigma_min_db, 'b', w, sigma_max_db, 'b')
    title ('Singular Values')
    xlabel ('Frequency [rad/s]')
    ylabel ('Singular Values [dB]')
    grid on
  else # Return values
    sigma_min_r = sigma_min;
    sigma_max_r = sigma_max;
    w_r = w;
  endif

endfunction


%!shared sigma_min_exp, sigma_max_exp, w_exp, sigma_min_obs, sigma_max_obs, w_obs
%! A = [1 2; 3 4];
%! B = [5 6; 7 8];
%! C = [4 3; 2 1];
%! D = [8 7; 6 5];
%! w = [2 3];
%! sigma_min_exp = [0.698526948925716   0.608629874340667];
%! sigma_max_exp = [7.91760889977901   8.62745836756994];
%! w_exp = [2 3];
%! [sigma_min_obs, sigma_max_obs, w_obs] = svplot (ss (A, B, C, D), w);
%!assert (sigma_min_obs, sigma_min_exp, 4*eps); # tolerance manually tweaked
%!assert (sigma_max_obs, sigma_max_exp, 12*eps); # tolerance manually tweaked
%!assert (w_obs, w_exp, 2*eps);