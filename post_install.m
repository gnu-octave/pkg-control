function post_install (d)

  fprintf ('\n');
  fprintf ('The control package was installed into the directory\n');
  fprintf ('%s.\n', d.dir);
  fprintf ('\n');
  fprintf ('This package is released under the GPLv3+, which can be\n');
  fprintf ('found in \"packinfo/COPYING\". License and copyright information\n');
  fprintf ('of the used SLICOT-Reference routines (BSD 3-Clause) are located\n');
  fprintf ('in \"doc/SLICOT\".\n');
  fprintf ('\n');
  fprintf ('The online documentation of the control package at\n');
  fprintf ('https://gnu-octave.github.io/pkg-control/ can be accessed\n');
  fprintf ('from within Octave with the command doc_control ().\n');
  fprintf ('\n');

end
