## Copyright (C) 2017  Mark Bronsfeld
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
## @deftypefn {Function File} damp(@var{sys})
## @deftypefnx {Function File} {[@var{Wn}, @var{zeta}] =} damp(@var{sys})
## @deftypefnx {Function File} {[@var{Wn}, @var{zeta}, @var{P}] =} damp(@var{sys})
## Calculate natural frequencies, damping ratios and poles.
##
## If no output is specified, display overview table containing poles, 
## magnitude (if @var{sys} is a discrete-time model), damping ratios, natural 
## frequencies and time constants.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} model.
## @end table
##
## @strong{Outputs}
## @table @var
## @item Wn
## Natural frequencies of each pole of @var{sys} (in increasing order).
## The frequency unit is rad/s (radians per second).
## If @var{sys} is a discrete-time model with specified sample time, @var{Wn} 
## contains the natural frequencies of the equivalent continuous-time poles 
## (see Algorithms).  If @var{sys} has an unspecified sample time 
## (@var{tsam =} -1), @var{tsam =} 1 is used to calculate @var{Wn}.
## @item zeta
## Damping ratios of each pole of @var{sys} (in the same order as @var{Wn}).
## If @var{sys} is a discrete-time model with specified sample time, 
## @var{zeta} contains the damping ratios of the equivalent continuous-time 
## poles (see Algorithms).  If @var{sys} has an unspecified sample time
## (@var{tsam =} -1), @var{tsam =} 1 is used to calculate @var{zeta}.
## @item P
## Poles of @var{sys} (in the same order as @var{Wn}).
## @end table
##
## @strong{Algorithm}@*
## @table @var
## @item Pole
## Poles s (or z for discrete-time models) are calculated via @command{pole} 
## and resorted in order of increasing natural frequency.
## @item Equivalent continuous-time pole
## s = log (z) / sys.tsam (discrete-time models only)
## @item Magnitude
## mag = abs (z) (discrete-time models only)
## @item Natural Frequency
## Wn = abs (s)
## @item Damping ratio
## zeta = -cos (arg (s))
## @item Time constant
## tau = 1 / (Wn * zeta)
## @end table
##
## @seealso{dsort, eig, esort, pole, pzmap, zero}
## @end deftypefn

## Author: Mark Bronsfeld <m.brnsfld@googlemail.com>
## Created: January 2017
## Version: 0.1

function [Wn_out, zeta, P] = damp (sys)

  if (nargin == 1) # damp (sys)
    if (! (isa (sys, "lti") || issquare (sys)))
      error ("damp: argument must be an LTI system");
    endif
  else
    print_usage ();
  endif
        
  P = pole (sys); # Poles
  
  ## Distinguish between system/state matrices, continuous- and ...
  ## discrete-time models
  if ((! (isa (sys, "lti")) && issquare (sys)) || (isct (sys)))
    s = P;
  elseif (isdt (sys))
    if (sys.tsam == -1) # If sample time is unspecified...
      s = log (P);      # ...assume 1 second: log (P) / 1
    else
      s = log (P) ./ sys.tsam;
    endif
    
    mag = abs (P);  # Magnitude
  endif
  
  Wn = abs (s); # Frequencies (rad / seconds)
  
  ## Sort all vectors in order of increasing natural frequency
  [Wn, ndx] = sort (Wn);
  P = P (ndx);
  s = s (ndx);
  
  zeta = -cos (arg (s));  # Damping
  
  tau = 1 ./ (Wn .* zeta);  # Time constant (seconds)
  
  ## Suppress "ans" output if no output specified (only assign "Wn_out" ...
  ## if any output specified)
  if (nargout > 0)
    Wn_out = Wn;
  ## Display overview table when no output specified
  elseif (nargout == 0)
    ## Type conversion and formatting to exponential format
    P =     num2str (P, '%1.2e');
    zeta =  num2str (zeta, '%1.2e');
    Wn =    num2str (Wn, '%1.2e');
    tau =   num2str (tau, '%1.2e');
    
    ## Construct columns of overview table
    Pole =          [['Pole'; ' '; ' ']; P];
    Damping =       [['Damping'; ' '; ' ']; zeta];
    Frequency =     [['Frequency'; '(rad/seconds)'; ' ']; Wn];
    TimeConstant =  [['Time Constant'; '(seconds)'; ' ']; tau];
    
    ## Construct overview table - distinguish between system/state ...
    ## matrices, continuous- and discrete-time models
    if ((! (isa (sys, "lti")) && issquare (sys)) || (isct (sys)))
      ## Overview table
      overview_table = [repmat(' ',rows(Pole),3) Pole ...
                        repmat(' ',rows(Pole),4) Damping ...
                        repmat(' ',rows(Pole),4) Frequency ...
                        repmat(' ',rows(Pole),4) TimeConstant];
    elseif (isdt (sys))
      mag = mag (ndx); # Sort vector in order of increasing natural frequency
      mag = num2str (mag, '%1.2e'); # Type conversion and formatting to ...
                                    # exponential format
      
      Magnitude = [['Magnitude'; ' '; ' ']; mag]; # Construct additional ...
                                                  # column of overview table
      
      ## Overview table
      overview_table = [repmat(' ',rows(Pole),3) Pole ...
                        repmat(' ',rows(Pole),4) Magnitude ...
                        repmat(' ',rows(Pole),4) Damping ....
                        repmat(' ',rows(Pole),4) Frequency ....
                        repmat(' ',rows(Pole),4) TimeConstant];
    endif
    
    disp (overview_table); # Display overview table
  endif

endfunction

## Test system/state matrix
%!test
%! A = [-1,  0,  0;
%!       0, -2,  0;
%!       0,  0, -3 ];
%!
%! Wn_exp = [1;
%!           2;
%!           3 ];
%!
%! zeta_exp = [1;
%!             1;
%!             1 ];
%!
%! P_exp = [-1;
%!          -2;
%!          -3 ];
%!
%! [Wn_obs, zeta_obs, P_obs] = damp (A);
%!
%! assert (Wn_obs, Wn_exp, 0);
%! assert (zeta_obs, zeta_exp, 0);
%! assert (P_obs, P_exp, 0);

## Test continuous-time model
%!test
%! H = tf ([2, 5, 1], [1, 2, 3]);
%!
%! Wn_exp = [1.7321;
%!           1.7321 ];
%!
%! zeta_exp = [0.5774;
%!             0.5774 ];
%!
%! P_exp = [-1.0000 + 1.4142i;
%!          -1.0000 - 1.4142i ];
%!
%! [Wn_obs, zeta_obs, P_obs] = damp (H);
%!
%! assert (Wn_obs, Wn_exp, 1e-4);
%! assert (zeta_obs, zeta_exp, 1e-4);
%! assert (P_obs, P_exp, 1e-4);

## Test discrete-time model
%!test
%! H = tf ([5, 3, 1], [1, 6, 4, 4], 0.01);
%!
%! Wn_exp = [193.4924;
%!           193.4924;
%!           356.5264 ];
%!
%! zeta_exp = [ 0.0774;
%!              0.0774;
%!             -0.4728 ];
%!
%! P_exp = [-0.3020 + 0.8063i;
%!          -0.3020 - 0.8063i;
%!          -5.3961 + 0.0000i ];
%!
%! [Wn_obs, zeta_obs, P_obs] = damp (H);
%!
%! assert (Wn_obs, Wn_exp, 1e-4);
%! assert (zeta_obs, zeta_exp, 1e-4);
%! assert (P_obs, P_exp, 1e-4);
