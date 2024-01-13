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
## Otherwise, the poles and zeros are computed and returned.
##
## @strong{Inputs}
## @table @var
## @item sys
## @acronym{LTI} model.
## @item 'style'
## Line style and color, e.g. 'r' for a solid red line or '-.k' for a dash-dotted
## black line.  See @command{help plot} for details.
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

  sys_idx = cellfun (@isa, varargin, {"lti"});      # look for LTI models
  sty_idx = cellfun (@ischar, varargin);            # look for strings (plot styles)

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
    pol_re = cellfun (@real, pol, "uniformoutput", false);
    pol_im = cellfun (@imag, pol, "uniformoutput", false);
    zer_re = cellfun (@real, zer, "uniformoutput", false);
    zer_im = cellfun (@imag, zer, "uniformoutput", false);

    ## extract plotting styles
    tmp = cumsum (sys_idx);
    tmp(sys_idx | ! sty_idx) = 0;
    n = nnz (sys_idx);
    sty = arrayfun (@(x) varargin(tmp == x), 1:n, "uniformoutput", false);

    colororder = get (gca, "colororder");
    rc = rows (colororder);
    def_pol = arrayfun (@(k) {"x", "linewidth", 2, "color", colororder(1+rem (k-1, rc), :)}, 1:n, "uniformoutput", false);
    def_zer = arrayfun (@(k) {"o", "linewidth", 2, "color", colororder(1+rem (k-1, rc), :)}, 1:n, "uniformoutput", false);
    idx_no_sty = cellfun (@isempty, sty);
    sty_pol = sty_zer = sty;
    sty_pol(idx_no_sty) = def_pol(idx_no_sty);
    sty_zer(idx_no_sty) = def_zer(idx_no_sty);

    leg_args = cell (1, n);
    idx = find (sys_idx);
    dt = false;
    for k = 1 : n
      if (! idx_no_sty(k))
        ## style given, only allow custom colors, no custom markers
        [opt,vopt] = __pltopt__ ('pzmap', sty{k}, false);
        if (! @isempty (opt.color))
          sty_pol{1,k} = {"x", "linewidth", 2, "color", opt.color};
          sty_zer{1,k} = {"o", "linewidth", 2, "color", opt.color};
        else
          if (! vopt)
            warning ("pzmap: ignoring undefined color value in style \"%s\"\n", sty{k}{1,1});
          endif
          sty_pol(k) = def_pol(k);
          sty_zer(k) = def_zer(k);
        endif
      endif
      dt = dt || (varargin{idx(k)}.tsam != 0);
    endfor

    pol_args = horzcat (cellfun (@horzcat, pol_re, pol_im, sty_pol, "uniformoutput", false){:});
    zer_args = horzcat (cellfun (@horzcat, zer_re, zer_im, sty_zer, "uniformoutput", false){:});

    hold on;

    ## If no zeroes then just plot the poles and vice versa
    h = [];
    leg_args = cell ();
    for k = 1:n
      name = inputname (idx(k));
      if (isempty (name))
        name = ["Sys ", num2str(k)];       # needed for  pzmap (lticell{:})
      endif
      if (isempty (zer_re{1,k}))
        hx = plot (pol_re{1,k}, pol_im{1,k}, sty_pol{1,k}{:});
        leg_args = { leg_args{:}, ["poles ", name] }; 
      elseif  (isempty (pol_re{1,k}))
        hx = plot (zer_re{1,k}, zer_im{1,k}, sty_zer{1,k}{:});
        leg_args = { leg_args{:}, ["zeros ", name] }; 
      else
        hx = plot (pol_re{1,k}, pol_im{1,k}, sty_pol{1,k}{:}, ...
                   zer_re{1,k}, zer_im{1,k}, sty_zer{1,k}{:});
        leg_args = { leg_args{:}, ["poles ", name], ["zeros ", name] }; 
      endif
      h = [ h; hx ];
    endfor

    grid ("on")  
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

    a = gca ();
    % avoid flickering while drawing axis / stablity region
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

    legend (h, leg_args)


    else
    pol_r = pol{1};
    zer_r = zer{1};
  endif
  
endfunction

%!demo
%! s = tf('s');
%! g = (s-1)/(s-2)/(s-3);
%! pzmap(g);

%!test
%! s = tf('s');
%! g = (s-1)/(s-2)/(s-3);
%! [pol zer] = pzmap(g);
%! assert(sort(pol), [2 3]', 2*eps);
%! assert(zer, 1, eps);

%!demo
%! s = tf('s');
%! g = 1/(2*s^2+3*s+4);
%! pzmap(g);

%!test
%! s = tf('s');
%! g = 1/(2*s^2+3*s+4);
%! pol = pzmap(g);
%! assert(sort(pol), sort(roots([2 3 4]')), eps);
