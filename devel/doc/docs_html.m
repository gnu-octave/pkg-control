## Copyright 2025 Torsten Lilge <ttl-octave@mailbox.org>
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

## Function for creating the html documentation for the Control Package

function docs_html (pkg_name)

  # load package and get it's installation dir
  pkg ("load", pkg_name);
  info = pkg ("list", pkg_name);
  inst_dir = info{1}.dir;

  # we are in devel/doc
  docs_html_dir = "../../docs";
  logo_dev = fullfile (pwd, [pkg_name, "-logo.svg"]);
  logo_doc = fullfile (inst_dir, "doc", [pkg_name, ".svg"]);

  # does logo file exist?
  if (! exist (logo_dev, "file"))
    error ("logo file %s does not exist", logo_dev);
  endif

  # copy logo file to doc dir of installed package
  # (ny this, it does not have to be in the tarball)
  st = copyfile (logo_dev, logo_doc);
  if (! st)
    error ("cannot copy %s to %s", logo_dev, logo_doc);
  endif

  # load pkg doc and create docs 
  pkg ("load", "pkg-octave-doc");
  cd (docs_html_dir);
  package_texi2html (pkg_name);

  # remove logo from installation dir
  delete (logo_doc);

endfunction
