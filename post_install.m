function post_install (d)

  fprintf ('\n');
  fprintf ('The control package was installed into the directory\n');
  fprintf ('%s.\n', d.dir);
  fprintf ('License information can be found in \"packinfo/COPYING\"\n');
  fprintf ('and for the used SLICOT-Reference library in \"doc/SLICOT\".\n');
  fprintf ('\n');
  fprintf ('The online documentation of the control package at\n');
  fprintf ('https://gnu-octave.github.io/pkg-control/ can be accessed\n');
  fprintf ('from within Octave with the command doc_control ().\n');
  fprintf ('\n');

end
