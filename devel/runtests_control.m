## This file could replace test_control.m 
## if the following bug gets fixed:
## http://savannah.gnu.org/bugs/?40053


dir = fileparts (which ("test_control"));

runtests (dir)
