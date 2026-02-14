## Copyright (C) 2023-2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2026      Torsten Lilge <ttl-octave@mailbox.org>
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.
##
## This code is taken from the package "pkg-octave-doc" and was adapted
## in order to add demo examples to the pdf and qch documentation of the
## control package

## -*- texinfo -*-
## @deftypefn  {Function File} {} build_demos (@var{pkgname}, @var{fcnname}, @var{texi_file})
##
## Generate texi code from the available demos of a particular function.
##
## ## @strong{Inputs}
## @table @var
## @item pkgname
## A char string with the package's name.
## @item fcnname
## A char string with the function's name.
## @item texi_file
## A char string with the name of the file where the resulting texinfo
## code of the demo has to be saved.
## @end table
## @end deftypefn

function build_demos (pkgname, fcnname, texi_file)

  if (nargin != 3)
    print_usage ();
  endif

  if (! ischar (fcnname))
    print_usage ();
  endif

  ## Get available demos from function
  pkg ("load", pkgname);
  texi = "";
  demos = {};

  [code, idx] = test (fcnname, "grabdemo");
  if (idx > 0)
    for i = 1:length (idx) - 1
      block = code(idx(i):idx(i+1)-1);
      demos = [demos; sprintf("%s\n", block)];
    endfor
  endif

  ## For @class methods: Clean up fileprefix
  fcnfile = strrep (fcnname, filesep, "_");

  if (! isempty (demos))

    ## Demos template
    demos_template = ["@strong{Example: {{NUMBER}}}\n", ...
                      "@example\n", ...
                      "@group\n", ...
                      "{{DEMO}}\n", ...
                      "@end group\n", ...
                      "@end example\n", ...
                      "{{FIGURES}}"];

    ## For each demo
    for demo_num = 1:numel (demos)
      try
        ## Initialize texi string for this demo
        demo_texi = "";

        ## Prepare environment variables
        close all
        diary_file = "__diary__.txt";
        if (exist (diary_file, "file"))
          delete (diary_file);
        endif

        unwind_protect

          ## Get current values
          dfv = get (0, "defaultfigurevisible");
          set (0, "defaultfigurevisible", "off");
          oldpager = PAGER('/dev/null');
          oldpso = page_screen_output(1);
          oldwarn = warning ();

          warning ("off");  # prevent any warning in demo output

          ## Format texi string with demo code
          demo_code = demos{demo_num};
          demo_texi = [demo_texi demo_code];

          ## Evaluate DEMO code
          diary (diary_file);
          eval (demo_code);
          diary ("off");

          ## Read __diary__.txt
          fid = fopen (diary_file);
          demo_text = fscanf (fid, "%c", Inf);
          fclose (fid);

          ## Format texi string with demo output
          if (! strcmp (demo_texi, "\n"))
            demo_texi = [demo_texi demo_text];
          endif

          ## Save figures
          images = {};
          figure_num = demo_num * 100;
          while (! isempty (get (0, "currentfigure")))
            figure_num = figure_num + 1;
            fig = gcf ();
            name = sprintf ("%s_%d", fcnfile, figure_num);
            fullpath = fullfile ("figures", name);
            ## use pdf for the pdf and png for the qch (html) documentation
            print (fig, [fullpath, ".pdf"], "-dpdfcrop", "-F:16", "-S560,420");
            print (fig, [fullpath, ".png"], "-dpng", "-F:16", "-S360,270");
            images{end+1} = fullpath;
            close (fig);
          endwhile

          ## Reverse image list, since we got them latest-first
          images = images (end:-1:1);

          ## Add reference to image (if applicable). Only the basename is
          ## used. makeinfo pdf output will use pdf and html output png files.
          demo_fig = "";
          if (! isempty (images))
            for i = 1:numel (images)
              demo_fig = [demo_fig sprintf("@center @image{%s,240px}\n", images{i})];
            endfor
          endif

          ## Append demo demo_texi to texi using the template
          full_demo_texi = strrep (demos_template, "{{NUMBER}}", ...
                              sprintf ("%d", demo_num));
          full_demo_texi = strrep (full_demo_texi, "{{DEMO}}", ...
                              demo_texi(2:end-1));
          full_demo_texi = strrep (full_demo_texi, "{{FIGURES}}", ...
                              sprintf ("%s", demo_fig));
          while (full_demo_texi(end) == "\n")
            full_demo_texi = full_demo_texi(1:end-1);
          endwhile
          full_demo_texi = [full_demo_texi "\n"];

        unwind_protect_cleanup

          delete (diary_file);
          set (0, "defaultfigurevisible", dfv);
          PAGER(oldpager);
          page_screen_output(oldpso);
          warning (oldwarn);

        end_unwind_protect

        texi = [texi full_demo_texi];

        ## write texi_file
        fid = fopen (texi_file, "W");
        fprintf (fid, "%s", texi);
        fclose (fid);

      catch
        printf ("Unable to process demo %d from %s:\n %s\n", ...
                demo_num, fcnname, lasterr);
      end_try_catch
    endfor
  endif

endfunction
