## Copyright (C) 2022 Torsten Lilge
##
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {} zgrid {@ }
## @deftypefnx {Function File} {} zgrid on
## @deftypefnx {Function File} {} zgrid off
## @deftypefnx {Function File} {} zgrid (@var{z}, @var{w})
## @deftypefnx {Function File} {} zgrid (@var{hax}, @dots{})
## Display an grid in the complex z-plane.
##
## Control the display of z-plane grid with :
## @itemize
## @item zeta lines corresponding to damping ratios and
## @item omega lines corresponding to undamped natural frequencies
## @end itemize
##
## The function state input may be either @qcode{"on"} or @qcode{"off"}
## for creating or removing the grid. If omitted, a new grid is created
## when it does not exist or the visibility of the current grid is toggled.
##
## The zgrid will automatically plot the grid lines at nice values or
## at constant values specified by two arguments :
##
## @example
## zgrid (@var{Z}, @var{W})
## @end example
##
## @noindent
## where @var{Z} and @var{W} are :
## @itemize
## @item @var{Z} vector of constant zeta values to plot as lines
## @item @var{W} vector of constant omega values to plot as circles
## @end itemize
##
## Example of usage:
## @example
## zgrid on	  create the z-plane grid
## zgrid off 	remove the z-plane grid
## zgrid		  toggle the z-plane grid visibility
## zgrid ([0.3, 0.8, @dots{}], [0.25*pid, 0.5*pi, @dots{}])   create:
## @example
## @itemize
## @item zeta lines for 0.3, 0.8, @dots{}
## @item omega lines for 0.25*pi/T, 0.5*pi/T, @dots{} [rad/s]
## @end itemize
## @end example
## zgrid (@var{hax}, @qcode{"on"})   create the z-plane grid for the axis
##                     handle @var{hax}
## @end example
##
## @seealso{grid,sgrid}
##
## @end deftypefn

## Author: Torsten Lilge based on "sgrid" by Stefan Mátéfi-Tempfli
## Created: 2023-05-14

function zgrid(varargin)

  [hax, varargin, nargs] = __plt_get_axis_arg__("zgrid", varargin{:});

  T = 1;

  if (nargs > 3)
    print_usage();
  endif

  if (isempty(hax))
    hax = gca();
  endif

  hg = findobj(hax, "tag", "zgrid");
  if (isempty(hg))
    hg = hggroup(hax, "tag", "zgrid");
    v_new = 1;
  else
        v_new = 0;
  endif

  v_z = []; v_w = [];

  if (nargs == 0)
    if (v_new)
      __zgrid_create__(hax, hg, v_z, v_w)
    else
      __zgrid_toggle__(hg)
    endif
  elseif (nargs == 1)
    arg1 = varargin{1};
    if (! ischar(arg1))
      error("zgrid: argument must be a string");
    endif
    if (strcmpi(arg1, "off"))
      __zgrid_delete__(hg);
    elseif (strcmpi(arg1, "on"))
      if (v_new)
        __zgrid_create__(hax, hg, v_z, v_w)
      else
        v_user = get(hg, "userdata");
        if (!isempty(v_user.z) | !isempty(v_user.w))
          __zgrid_delete_handles__(hg);
          __zgrid_create__(hax, hg, v_z, v_w)
        elseif (strcmp(get(hg, "visible"), "off"))
          set(hg, "visible", "on");
        endif
      endif
    else
      print_usage();
    endif
  else
    v_z = varargin{1};
    if (! isnumeric(v_z))
      error ("zgrid: Z argument (1) must be numeric");
    endif
    if (any(v_z < 0 | v_z > 1))
      error("zgrid: Z argument (1) must have values betwenn 0 .. 1");
    endif
    v_w = varargin{2};
    if (! isnumeric(v_w))
      error("zgrid: W argument (2) must be numeric");
    endif
    if (any(v_w < 0))
      error("zgrid: W argument (2) must have values larger or equal 0");
    endif
    if (any(v_w > pi/T))
      error("zgrid: W argument (2) must have values smaller or equal pi");
    endif
    if (v_new)
      __zgrid_create__(hax, hg,v_z, v_w)
    else
      v_user = get(hg, "userdata");
      if (!isequal(v_z, v_user.z) || !isequal(v_w, v_user.w) || (v_user.cl != [xlim() ylim()]))
        __zgrid_delete_handles__(hg);
        __zgrid_create__(hax, hg, v_z, v_w)
      elseif (strcmp(get(hg, "visible"), "off"))
        set(hg, "visible", "on");
      endif
    endif
  endif
endfunction

##----------------------------------------------------
function __zgrid_create__(hax, hg, v_z, v_w)

    hold on;
    box on;
    v_user.z = v_z;
    v_user.w = v_w;
    v_user.cl = [ xlim() ylim() ];

    set(hg, "userdata", v_user);

    T = 1;
    w_max = pi/T;
    d_w = pi/T/10;
    if (isempty(v_w))
      v_w = d_w:d_w:w_max;
    endif

    d_D = (v_user.cl(4)-v_user.cl(3))/50;
    D_max = 1;
    D = 0:d_D:D_max;
    clear ('i');

    for iw = 1:length(v_w)

      z = exp (T*(-D.*v_w(iw))) .* exp (T*i*v_w(iw).*sqrt (1-D.^2));
      z = [ conj(z) flip(z) ];
      plot (z, ":k", "linewidth", 0.6, "parent", hg);

      z_vis = __zgrid_visible__(z,v_user.cl);
      if (! isempty (z_vis))
        num = sprintf ("%1.1f \\pi/T", v_w(iw)*T/pi);
        zt = z_vis(end); % exp (T*i*v_w(iw));
        text (real (zt), imag (zt), num, "parent", hg);
      endif

    endfor

    d_D = 0.1;
    d_w = pi/T/50;
    if (isempty(v_z))
      v_z = 0:d_D:D_max-eps;
    endif

    clear ('i');
    wt = 0.25*pi/T;

    for (iz = 1:length(v_z))

      % Compute w for reaching the negative real axis (pi/T for D = 1):
      % T*w_end*sqrt(1-D^2) == pi
      w_end = pi/T/sqrt(1-v_z(iz)^2);
      w = [ 0:d_w:w_max w_max:(w_end-w_max)/10:w_end ];
      z = exp (T*(-v_z(iz).*w)) .* exp (i*T*w.*sqrt (1-v_z(iz).^2));

      plot (z,":k", "linewidth", 0.6, "parent", hg);
      plot (conj (z), ":k", "linewidth", 0.6, "parent", hg);

      z_vis = __zgrid_visible__(z,v_user.cl);
      if (! isempty (z_vis))
        num = sprintf ("%1.1f", v_z(iz));
        zt_idx = find (imag (z_vis) == max (imag (z_vis)));
        if zt_idx > 3
          zt_idx = zt_idx - 3;
        else
          zt_idx = 1;
        endif
        zt = z_vis (zt_idx);
        text (real (zt), imag (zt), num, "parent", hg);
      endif

    endfor

    hold off;
endfunction

##----------------------------------------------------
function __zgrid_delete__(hg)
  delete(hg);
endfunction

##----------------------------------------------------
function __zgrid_delete_handles__(hg)
  delete(get(hg, "children"));
endfunction

##----------------------------------------------------
function __zgrid_toggle__(hg)
  if (strcmp(get(hg, "visible"), "on"))
    set(hg, "visible", "off");
  elseif (strcmp(get(hg, "visible"), "off"))
    set(hg, "visible", "on");
  endif
endfunction

function z_vis = __zgrid_visible__ (z, cl)
  z_vis = z;
  z_vis = z_vis(real (z_vis) > cl(1));
  z_vis = z_vis(real (z_vis) < cl(2));
  z_vis = z_vis(imag (z_vis) > cl(3));
  z_vis = z_vis(imag (z_vis) < cl(4));
endfunction

%!demo
%! clf;
%! num = [1 0.25]; den = [1 -1.5 0];
%! sys = tf(num, den, 1);
%! rlocus(sys,0.01,0,3.5);
%! ylim([-1.1,1.1]);
%! zgrid on;

