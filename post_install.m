function post_install (d)

  fprintf ('\n');
  fprintf ('The control package was installed into the directory\n');
  fprintf ('%s.\n', d.dir);
  fprintf ('License and copyright information can be found in ');
  fprintf ('\"packinfo/COPYING\".\n');

  fprintf ('\n');
  fprintf ('If the control package was updated in this Octave\n');
  fprintf ('session, you might have to restart Octave.\n');
  fprintf ('\n');

end
