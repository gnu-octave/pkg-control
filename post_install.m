function post_install (d)

  fprintf ('\n');
  fprintf ('The control package was installed into the directory <ID>:\n');
  fprintf ('%s.\n', d.dir);
  fprintf ('  - License and copyright: \"<ID>/packinfo/COPYING\"\n');
  fprintf ('  - Documentation:\n');
  fprintf ('    - PDF:    \"<ID>/doc/control.pdf\"\n');
  fprintf ('    - QtHelp: Octave GUI documentation browser\n');
  fprintf ('    - Online: https://gnu-octave.github.io/pkg-control\n');
  fprintf ('\n');
  fprintf ('If the control package was updated in this Octave\n');
  fprintf ('session, you might have to restart Octave.\n');
  fprintf ('\n');

end
