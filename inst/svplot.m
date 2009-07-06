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
## Version: 0.1alpha


function [sigma_min_r, sigma_max_r, w_r] = svplot (sys, w)

  if (nargin < 1 || nargin > 2)
    print_usage();
  endif

  ## Get State Space Matrices
  sys = sysupdate(sys, "ss");
  [A, B, C, D] = sys2ss(sys); 


  ## Find interesting Frequency Range w if not specified
  if (nargin < 2)

    j = 1; # Index Counter
    
    for m = 1 : size(B, 2) # For all Inputs
      for p = 1 : size(C, 1) # For all Outputs

        ## sys2zp and bode_bounds don't work for MIMO Systems!
        sysp = sysprune(sys, p, m);
        [zer, pol, k, tsam, inname, outname] = sys2zp(sysp);

        if (tsam == 0)
          DIGITAL = 0;
        else
          DIGITAL = 1;
        endif

        [wmin, wmax] = bode_bounds(zer, pol, DIGITAL, tsam);

        w_min(j) = wmin;
        w_max(j) = wmax;
        j++;
      endfor
    endfor

    dec_min = min(w_min); # Begin Plot at 10^dec_min rad/s
    dec_max = max(w_max); # End Plot at 10^dec_min rad/s
    n_steps = 1000; # Number of Steps
    
    w = logspace(dec_min, dec_max, n_steps); # rad/s
  endif


  ## Repeat for every Frequency w
  for k = 1 : size(w, 2)

    ## Frequency Response Matrix
    P = C * inv(i * w(k) * eye(size(A)) - A) * B  +  D;

    ## Singular Value Decomposition
    sigma = svd(P);
    sigma_min(k) = min(sigma);
    sigma_max(k) = max(sigma);
  endfor

  if (nargout < 1) # Plot the Information

    ## Convert to dB for Plotting
    sigma_min_db = 20 * log10(sigma_min);
    sigma_max_db = 20 * log10(sigma_max);
   
    ## Plot Results
    semilogx(w, sigma_min_db, 'b', w, sigma_max_db, 'b')
    title('Singular Values')
    xlabel('Frequency [rad/s]')
    ylabel('Singular Values [dB]')
    grid on
  else # Return Values
    sigma_min_r = sigma_min;
    sigma_max_r = sigma_max;
    w_r = w;
  endif

endfunction
