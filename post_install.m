function post_install (d)

  fprintf ('The control package uses some routines of the ');
  fprintf ('SLICOT-Reference library (https://github.com/SLICOT/SLICOT-Reference), ');
  fprintf ('which are available under the BSD 3-Clause License, ');
  fprintf ('which can be found in %s/doc/slicot together with a ', d.dir);
  fprintf ('readme file.\n');

end
