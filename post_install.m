function post_install (d)

  fprintf ('\n');
  fprintf ('The control package was installed into the directory\n');
  fprintf ('%s.\n', d.dir);
  fprintf ('License and copyright information can be found in ');
  fprintf ('\"packinfo/COPYING\"\n');

  fprintf ('\n');
  fprintf ('Please report bugs and feature requests at\n');
  fprintf ('https://github.com/gnu-octave/pkg-control/issues\n');

  fprintf ('\n');
  fprintf ('If the control package was updated in the same Octave\n');
  fprintf ('session, you might have to restart Octave.\n');

  fprintf ('\n');

end
