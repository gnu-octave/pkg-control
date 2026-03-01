## Copyright (C) 2009-2016   Lukas F. Reichlin
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
## @deftypefn {Function File} {} pzmap (@var{sys})
## @deftypefnx {Function File} {} pzmap (@var{sys1}, @var{sys2}, @dots{}, @var{sysN})
## @deftypefnx {Function File} {} pzmap (@var{sys1}, @var{'style1'}, @dots{}, @var{sysN}, @var{'styleN'})
## @deftypefnx {Function File} {[@var{p}, @var{z}] =} pzmap (@var{sys})
## Plot the poles and zeros of an LTI system in the complex plane.
## If no output arguments are given, the result is plotted on the screen.
## Otherwise, the poles and zeros are computed and returned. Note that
## only one system is processed when output arguments are given.
##
## @strong{Inputs}
## @table @var
## @item sys, sys1, ...
## @acronym{LTI} model(s).
## @item 'style'
## Color, e.g. 'r' for a red. See @command{help plot} for details.
## Marker or line styles are ignored as poles and zeros have the
## fixed marker types 'x' and 'o' repectively.
## @end table
##
## @strong{Outputs}
## @table @var
## @item p
## Poles of @var{sys}.
## @item z
## Invariant zeros of @var{sys}.
## @end table
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: November 2009
## Version: 0.2

function [pol_r, zer_r] = pzmap (varargin)

  if (nargin == 0)
    print_usage ();
  endif

  names = arrayfun (@inputname, 1:nargin, 'uniformoutput', false);

  [sys_idx, sty_idx, names, sty] = __control_args__ (varargin, {"@lti", @ischar}, names);

  inv_idx = ! (sys_idx | sty_idx);                  # invalid arguments

  if (any (inv_idx))
    warning ("pzmap: arguments number %s are invalid and are being ignored\n", ...
             mat2str (find (inv_idx)(:).'));
  endif

  if (nnz (sys_idx) == 0)
    error ("pzmap: require at least one LTI model");
  endif

  if (nargout > 0 && (nnz (sys_idx) > 1 || any (sty_idx)))
    print_usage ();
  endif

  if (any (find (sty_idx) < find (sys_idx)(1)))
    warning ("pzmap: strings in front of first LTI model are being ignored\n");
  endif

  pol = cellfun (@pole, varargin(sys_idx), "uniformoutput", false);
  zer = cellfun (@zero, varargin(sys_idx), "uniformoutput", false);

  if (! nargout)

    ms = 10;
    if (strcmp(graphics_toolkit(),'gnuplot'))
      ms = 6;
    endif

    pol_re = cellfun (@real, pol, "uniformoutput", false);
    pol_im = cellfun (@imag, pol, "uniformoutput", false);
    zer_re = cellfun (@real, zer, "uniformoutput", false);
    zer_im = cellfun (@imag, zer, "uniformoutput", false);

    ## extract plotting styles
    sty_pol = sty_zer = sty;
    idx = find (sys_idx);
    dt = false;
    for k = 1 : length(sty)
      style = find (cellfun (@(s) ischar(s) && (! strcmp (s, "color")), sty{k}));
      col = find (cellfun (@is_matrix, sty{k}));
      [opt,vopt] = __pltopt__ ('pzmap', sty{k}{style}, false);
      if (isempty (opt.color))
        c = sty{k}{col};
      else
        c = opt.color;
      endif
      sty_pol{1,k}{style} = ["x;poles ",opt.key,";"];
      sty_zer{1,k}{style} = ["o;zeros ",opt.key,";"];
      sty_pol{1,k} = {sty_pol{1,k}{:}, "linewidth", 2, "markersize", ms, "color", c};
      sty_zer{1,k} = {sty_zer{1,k}{:}, "linewidth", 2, "markersize", ms, "color", c};
      dt = dt || (varargin{idx(k)}.tsam != 0);
    endfor

    pol_args = horzcat (cellfun (@horzcat, pol_re, pol_im, sty_pol, "uniformoutput", false){:});
    zer_args = horzcat (cellfun (@horzcat, zer_re, zer_im, sty_zer, "uniformoutput", false){:});

    hold on;

    ## If no zeroes then just plot the poles and vice versa
    h = [];
    leg_args = cell ();
    for k = 1:length (sty)
      if (isempty (zer_re{1,k}))
        hx = plot (pol_re{1,k}, pol_im{1,k}, sty_pol{1,k}{:});
      elseif  (isempty (pol_re{1,k}))
        hx = plot (zer_re{1,k}, zer_im{1,k}, sty_zer{1,k}{:});
      else
        hx = plot (pol_re{1,k}, pol_im{1,k}, sty_pol{1,k}{:}, ...
                   zer_re{1,k}, zer_im{1,k}, sty_zer{1,k}{:});
      endif
      h = [ h; hx ];
    endfor

    legend ("autoupdate", "off");

    grid on
    box on
    title ("Pole-Zero Map")
    xlabel ("Real Axis")
    ylabel ("Imaginary Axis")

    xl = xlim();
    yl = ylim();

    dx = (xl(2)-xl(1))/10;
    xl(1) = xl(1) - dx;
    xl(2) = xl(2) + dx;
    dy = (yl(2)-yl(1))/10;
    yl(1) = yl(1) - dy;
    yl(2) = yl(2) + dy;

    ## avoid flickering while drawing axis / stablity region
    a = gca ();
    set (a, 'xlimmode', 'manual');
    set (a, 'ylimmode', 'manual');
    plot ([0,0], [yl(1)-dy*100, yl(2)+dy*100], '-', 'color', [0.7 0.7 0.7]);
    plot ([xl(1)-dx*100, xl(2)+dx*100], [0,0], '-', 'color', [0.7 0.7 0.7]);
    if dt
      t = 0:0.05:6.3;
      plot(cos(t),sin(t),'-', 'color', [0.7 0.7 0.7]);
    endif
    set (a, 'xlimmode', 'auto');
    set (a, 'ylimmode', 'auto');

    xlim(xl);
    ylim(yl);

    hold off;

  else

    ## output arguments: only one system is allowed as input
    pol_r = pol{1};
    zer_r = zer{1};

  endif

endfunction


%!demo
%! z = tf('z',1);
%! G1z = (z+1)/(z-0.75)/(z^2-1*z+1);
%! pzmap(G1z);

%!demo
%! s = tf('s');
%! G1 = 1/(2*s^2+3*s+4);
%! G2 = (1-s)/(1+s)/(s^2+s+1);
%! pzmap(G1,G2);

%!test
%! s = tf('s');
%! G = (s-1)/(s-2)/(s-3);
%! [pol zer] = pzmap(G);
%! assert(sort(pol), [2 3]', 2*eps);
%! assert(zer, 1, eps);

%!test
%! s = tf('s');
%! g = 1/(2*s^2+3*s+4);
%! pol = pzmap(g);
%! assert(sort(pol), sort(roots([2 3 4]')), eps);
