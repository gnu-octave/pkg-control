## Copyright (C) 2009, 2011   Lukas F. Reichlin
##
## This file is part of LTI Syncope.
##
## LTI Syncope is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## LTI Syncope is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with LTI Syncope.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} pzmap (@var{sys})
## @deftypefnx {Function File} {[@var{p}, @var{z}] =} pzmap (@var{sys})
## Plot the poles and zeros of an LTI system in the complex plane.
## If no output arguments are given, the result is plotted on the screen.
## Otherwise, the poles and zeros are computed and returned.
##
## @strong{Inputs}
## @table @var
## @item sys
## LTI model.
## @end table
##
## @strong{Outputs}
## @table @var
## @item p
## Poles of @var{sys}.
## @item z
## Transmission zeros of @var{sys}.
## @end table
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: November 2009
## Version: 0.1

function [pol_r, zer_r] = pzmap2 (varargin)

  ## TODO: multiplot feature:   pzmap (sys1, "b", sys2, "r", ...)

  if (nargin == 0)
    print_usage ();
  endif

  % if (! isa (sys, "lti"))
  %   error ("pzmap: argument must be an LTI system");
  % endif

  pol = cellfun ("@lti/pole", varargin, "uniformoutput", false);
  zer = cellfun ("@lti/zero", varargin, "uniformoutput", false);

  if (! nargout)
    pol_re = cellfun (@real, pol, "uniformoutput", false);
    pol_im = cellfun (@imag, pol, "uniformoutput", false);
    zer_re = cellfun (@real, zer, "uniformoutput", false);
    zer_im = cellfun (@imag, zer, "uniformoutput", false);
    
    len = numel (pol);
    plot_args = {};
    legend_args = cell (len, 1);
    colororder = get (gca, "colororder");
    for k = 1 : len
      plot_args = cat (2, plot_args, pol_re(k), pol_im(k), {"x", "color", colororder(k,:)}, \
                                     zer_re(k), zer_im(k), {"o", "color", colororder(k,:)});
      legend_args{k} = inputname(k);
    endfor

    h = plot (plot_args{:});
    grid ("on")  
    title ("Pole-Zero Map")
    xlabel ("Real Axis")
    ylabel ("Imaginary Axis")
    legend (h(1:2:2*len), legend_args)
  else
    pol_r = pol{1};
    zer_r = zer{1};
  endif
  
endfunction
