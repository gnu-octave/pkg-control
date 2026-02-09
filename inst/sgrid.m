## Copyright (C) 2019 Stefan Mátéfi-Tempfli
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
## @deftypefn  {Function File} {} sgrid {@ }
## @deftypefnx {Function File} {} sgrid on
## @deftypefnx {Function File} {} sgrid off
## @deftypefnx {Function File} {} sgrid (@var{z}, @var{w})
## @deftypefnx {Function File} {} sgrid (@var{hax}, @dots{})
##
## Display an grid in the complex s-plane.
##
## Control the display of s-plane grid with :
## @itemize
## @item zeta lines corresponding to damping ratios and
## @item omega circles corresponding to undamped natural frequencies
## @end itemize
##
## The function state input may be either @qcode{"on"} or @qcode{"off"}
## for creating or removing the grid. If omitted, a new grid is created
## when it does not exist or the visibility of the current grid is toggled.
##
## The sgrid will automatically plot the grid lines at nice values or
## at constant values specified by two arguments :
##
## @example
## sgrid (@var{Z}, @var{W})
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
## sgrid on	create the s-plane grid
## sgrid off 	remove the s-plane grid
## sgrid		toggle the s-plane grid visibility
## sgrid ([0.3, 0.8, @dots{}], [10, 75, @dots{}])   create:
##               zeta lines for 0.3, 0.8, @dots{}
##               omega circles for 10, 75, @dots{} [rad/s]
## sgrid off; sgrid  remove current s-grid and
##                   draw new with default values
## sgrid (@var{hax}, @qcode{"on"})   create the s-plane grid for the axis
##                     handle @var{hax}
## @end example
##
## @seealso{grid,zgrid}
##
## @end deftypefn

## Author: Stefan Mátéfi-Tempfli
## Created: 2019-01-24

function sgrid(varargin)

  [hax, varargin, nargs] = __plt_get_axis_arg__("sgrid", varargin{:});

  if (nargs > 3)
    print_usage();
  endif

  if (isempty(hax))
    hax = gca();
  endif

  hg = findobj(hax, "tag", "sgrid");
  if (isempty(hg))
    hg = hggroup(hax, "tag", "sgrid");
    v_new = 1;
  else
        v_new = 0;
  endif

  v_z = []; v_w = [];

  if (nargs == 0)
    if (v_new)
      __sgrid_create__(hax, hg, v_z, v_w)
    else
      __sgrid_toggle__(hg)
    endif
  elseif (nargs == 1)
    arg1 = varargin{1};
    if (! ischar(arg1))
      error("sgrid: argument must be a string");
    endif
    if (strcmpi(arg1, "off"))
      __sgrid_delete__(hg);
    elseif (strcmpi(arg1, "on"))
      if (v_new)
        __sgrid_create__(hax, hg, v_z, v_w)
      else
        v_user = get(hg, "userdata");
        if (!isempty(v_user.z) | !isempty(v_user.w))
          __sgrid_delete_handles__(hg);
          __sgrid_create__(hax, hg, v_z, v_w)
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
      error ("sgrid: Z argument (1) must be numeric");
    endif
    if (any(v_z < 0 | v_z > 1))
      error("sgrid: Z argument (1) must have values betwenn 0 .. 1");
    endif
    v_w = varargin{2};
    if (! isnumeric(v_w))
      error("sgrid: W argument (2) must be numeric");
    endif
    if (any(v_w <= 0))
      error("sgrid: W argument (1) must have positive values larger than 0");
    endif
    if (v_new)
      __sgrid_create__(hax, hg,v_z, v_w)
    else
      v_user = get(hg, "userdata");
      if (!isequal(v_z, v_user.z) || !isequal(v_w, v_user.w))
        __sgrid_delete_handles__(hg);
        __sgrid_create__(hax, hg, v_z, v_w)
      elseif (strcmp(get(hg, "visible"), "off"))
        set(hg, "visible", "on");
      endif
    endif
  endif
endfunction

##----------------------------------------------------
function __sgrid_create__(hax, hg, v_z, v_w)

    if nargin < 1 || isempty(hax)
    hax = gca();
  end

  hf = ancestor(hax, "figure");

  % Store handles for use in resize callback
  data.hax = hax;
  data.grid_state = "on";
  data.z = v_z;
  data.w = v_w;
  set(hg, "userdata", data);


  % Set resize callback
  set(hf, "resizefcn", @(src, evt) sgrid_resize_callback(src));


    hold on;
    box on;
    v_axis = axis(hax);
    if ((v_axis(2) < 0.15 * (v_axis(2) - v_axis(1))))
      v_axis(2) = 0.15 * (v_axis(2) - v_axis(1));
      axis(hax,v_axis);
    endif
    v_daxis = [v_axis(2) - v_axis(1), v_axis(4) - v_axis(3)];
    v_1d5axis = 0.015 * v_daxis;
    v_taxis = [v_axis(1) + v_1d5axis(1);
              v_1d5axis(1);
              v_axis(3) + v_1d5axis(2);
              v_axis(4) - v_1d5axis(2)];
    v_max = max (abs(v_axis)) * sqrt(2);
    v_rd = v_max/14;
    v_x = ceil(log10 (v_rd) - 1);
    v_p10 = 10^v_x;
    v_d = ceil(v_rd/v_p10) * v_p10;

    if (isempty(v_w))
      v_w = v_d:v_d:v_max;
    endif
    v_w_x = v_taxis(2) * ones(size(v_w));
    if ((0.1 < v_d) && (v_d < 1000))
      if (v_d >= 1)
        v_w_str = cellstr(int2str(v_w(:)));
      else
        v_w_str = cellstr(num2str(v_w(:),"%2.2f"));
      endif
    else
      v_w_str = cellstr(num2str(v_w(:),"%1.2e"));
    endif
    for (i = 1:length(v_w))
      v_sgrid.w(i) = plot(v_w(i)*cos(pi/2:0.01:3*pi/2),
                          v_w(i)*sin(pi/2:0.01:3*pi/2),
                          ":k", "linewidth", 0.6,
                          "parent", hg);
    endfor

    idx = sum(v_w <= v_taxis(4));
    for (i = 1:idx)
      text(v_w_x(i), v_w(i), v_w_str{i}, "parent", hg);
    endfor
    text(v_taxis(2), 0, "0 rad/s", "parent", hg);

    idx = sum(v_w <= -v_taxis(3));
    for (i = 1:idx)
      text(v_w_x(i), -v_w(i), v_w_str{i}, "parent", hg);
    endfor
    if (isempty(v_z))
      v_z = cos(linspace(0, pi/2, 13));   # every 7.5°
    endif
    v_z_phi = pi - [(acos(v_z(v_z != 0))), ...
                    (-acos(v_z((v_z != 0)&(v_z != 1))))];
    v_z_px = v_max * cos(v_z_phi);
    v_z_py = v_max * sin(v_z_phi);
    v_z_str = cellstr(num2str(abs(sin(v_z_phi(:) - pi/2)),"%1.2f"));
    for (i = 1:length(v_z_phi))
      plot([0 v_z_px(i)], [0 v_z_py(i)],
           ":k", "linewidth", 0.6,
           "parent", hg);
      pry = v_taxis(1) * tan(v_z_phi(i));
      if (pry > v_taxis(4))
        v_z_lc = v_taxis(4) * [cot(v_z_phi(i)) 1];
      elseif (pry < v_taxis(3))
        v_z_lc = v_taxis(3) * [cot(v_z_phi(i)) 1];
      else
        v_z_lc = v_taxis(1) * [1 tan(v_z_phi(i))];
      endif
      text (v_z_lc(1), v_z_lc(2), v_z_str{i}, "parent", hg);
    endfor
    plot([0 0], [0 v_max], "-k", "linewidth", 0.3, "parent", hg);
    plot([0 0], [0 -v_max], "-k", "linewidth", 0.3, "parent", hg);

    hold off;
endfunction

function sgrid_resize_callback(hf)
  % Get the grid handle (you might be storing it as data.hg or similar)
  hg = findobj(hf, "tag", "sgrid");  % Use your actual tag or method of storing the grid handle

  if isempty(hg)
    return;
  endif

  data = get(hg, "userdata");

  if isfield(data, "hax") && isgraphics(data.hax) && ...
   isfield(data, "v_z") && isfield(data, "v_w")
  __sgrid_create__(data.hax, hg, data.v_z, data.v_w);
  endif
endfunction


##----------------------------------------------------
function __sgrid_delete__(hg)
  delete(hg);
endfunction

##----------------------------------------------------
function __sgrid_delete_handles__(hg)
  delete(get(hg, "children"));
endfunction

##----------------------------------------------------
function __sgrid_toggle__(hg)
  if (strcmp(get(hg, "visible"), "on"))
    set(hg, "visible", "off");
  elseif (strcmp(get(hg, "visible"), "off"))
    set(hg, "visible", "on");
  endif
endfunction

%!demo
%! clf;
%! num = 1; den = [1 4 5 0];
%! sys = tf(num, den);
%! rlocus(sys);
%! sgrid on;

%!demo
%! clf;
%! num = [1 3]; den = [1 5 20 16 0];
%! sys = tf(num, den);
%! rlocus(sys);
%! hfig = get(0, "currentfigure");
%! hax = get(hfig, "currentaxes");
%! sgrid(hax, "on");

%!demo
%! clf;
%! num = 1; den = [1 4 5 0];
%! sys = tf(num, den);
%! rlocus(sys);
%! sgrid([0.5 0.7 0.89], [1 1.66 2.23]);

