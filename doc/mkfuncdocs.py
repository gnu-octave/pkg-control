#!/usr/bin/python3

## Copyright 2018-2024 John Donoghue
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

## mkfuncdocs v1.0.8
## mkfuncdocs.py will attempt to extract the help texts from functions in src
## dirs, extracting only those that are in the specifed INDEX file and output them
## to stdout in texi format
##
## It will extract from both .m and the help text for DEFUN_DLD help in .cc/.cpp
## files.
##
## It attempts to find the help text for each function in a file within the src search
## folders that match in order: [ functionname.m functionname.cc functionname.cpp
## functionname_withoutprefix.cc functionname_withoutprefix.cpp ]
##
## Usage:
##   mkfundocs.py options INDEXfile
## Options can be 0 or more of:
##   --verbose       : Turn on verbose mode
##   --src-dir=xxxxx : Add dir xxxxx to the dirs searched for the function file.
##                     If no directories are provided, it will default to looking in the
##                     'inst' directory.
##   --ignore=xxxxx  : dont attempt to generate help for function xxxxx.
##   --funcprefix=xxxxx : remove xxxxx from the function name when searching for matching
##                     source file.
##   --allowscan     : if can not find function, attemp to scan .cc,cpp,cxx files for match
##
##   --standalone    : generate a texinfo file expected to be used with being included in
##                     another document file.

import sys
import os
import re
import tempfile
import shutil
import fnmatch
import subprocess
import glob
import calendar;
import time;

class Group:
  name = "Functions"
  functions = []

  def __init__ (self, name=""):
    if name:
        self.name = name
    self.functions = []

class Index:
  name = ""
  groups = []

def texify_line(line):
  # convert any special chars in a line to texinfo format
  # currently just used for group formatting ?
  line = line.replace("@", "@@")
  line = line.replace("{", "@{")
  line = line.replace("}", "@}")
  line = line.replace(",", "@comma{}")
  return line

def find_defun_line_in_file(filename, fnname):
  linecnt = 0

  defun_line=re.compile(r"^DEFUN_DLD\s*\(\s*{}".format(fnname))
  with open(filename, 'rt') as f:
    for line in f:
      if re.match(defun_line, line):
        return linecnt

      linecnt = linecnt + 1

  return -1

def find_function_line_in_file(filename, fnname):
  linecnt = 0
  func = False
  defun_line=re.compile(r"^\s*function \s*")
  with open(filename, 'rt') as f:
    for line in f:
      if func == True:
          x = line.strip()
          if x.startswith("## -*- texinfo -*-"):
            return linecnt
          else:
            func = False

      if re.match(defun_line, line):
        if line.find("=") != -1:
           x = line.split("=")
           x = x[-1]
        else:
           x = line.replace("function ", "")

        x = x.split("(")
        x = x[0].strip()
        if x == fnname:
          func = True

      linecnt = linecnt + 1

  return -1


def read_m_file(filename, skip=0):
  help = []
  inhelp = False
  havehelp = False;
  with open(filename, 'rt') as f:
    for line in f:
      line = line.lstrip()
      if skip > 0:
        skip = skip - 1
      elif not havehelp:
        if havehelp == False and inhelp == False and line.startswith('##'):
          if "texinfo" in line:
            inhelp = True
        elif inhelp == True:
          if  not line.startswith('##'):
            inhelp = False
            havehelp = True
          else:
            if line.startswith("## @"):
              line = line[3:]
            else:
              line = line[2:]
            help.append (line.rstrip());

  return help

def read_cc_file(filename, skip=0):
  help = []
  inhelp = False
  havehelp = False;
  with open(filename, 'rt') as f:
    for line in f:
      line = line.lstrip()
      if skip > 0:
        skip = skip - 1
      elif not havehelp:
        if havehelp == False and inhelp == False:
          if "texinfo" in line:
            inhelp = True
        elif inhelp == True:
          line = line.rstrip()
          if len(line) > 0 and line[-1] == '\\':
            line = line[:-1]
            line = line.rstrip()

          line = line.replace("\\n", "\n") 
          line = line.replace("\\\"", "\"") 

          if len(line) > 0 and line[-1] == '\n':
            line = line[:-1]
          # endif a texinfo line
          elif line.endswith('")'):
            line = line[:-2]
 
          if  line.startswith('{'):
            inhelp = False
            havehelp = True
          else:
            help.append (line);

  return help

def read_help (filename, skip=0):
  help = []

  if filename[-2:] == ".m":
    help = read_m_file(filename, skip)
  else:
    help = read_cc_file(filename, skip)

  return help

def read_index (filename, ignore):
  index = Index ()

  with open(filename, 'rt') as f:
    lines = f.read().splitlines()

  #print ("read", lines)
  first = True
  category = Group()
  for l in lines:
    if l.startswith("#"):
      pass
    elif first:
      index.name = l;
      first = False
    elif l.startswith(" "):
        l = l.strip()
        # may be multiple functions here
        funcs = l.split()
        for f in funcs:
          if f not in ignore:
            category.functions.append(f);
    else:
      # new category name
      if len(category.functions) > 0:
        index.groups.append(category)
      category = Group(l.strip())

  # left over category ?
  if len(category.functions) > 0:
    index.groups.append(category)

  return index;

def find_class_file(fname, paths):
  
  for f in paths:
      # class constructor ?
      name = f + "/@" + fname + "/" + fname + ".m"
      if os.path.isfile(name):
        return name, 0

      # perhaps classname.func format ?
      x = fname.split(".")
      if len(x) > 0:
          zname = x.pop()
          cname = ".".join(x)
          name = f + "/" + cname + ".m"
          if os.path.isfile(name):
            idx = find_function_line_in_file(name, zname)
            if idx >= 0:
              return name, idx
          name = f + "/@" + cname + "/" + zname + ".m"
          if os.path.isfile(name):
            return name, 0
  return None, -1

def find_func_file(fname, paths, prefix, scanfiles=False):
  for f in paths:
      name = f + "/" + fname + ".m"
      if os.path.isfile(name):
        return name, 0
      # class constructor ?
      name = f + "/@" + fname + "/" + fname + ".m"
      if os.path.isfile(name):
        return name, 0
      name = f + "/" + fname + ".cc"
      if os.path.isfile(name):
        return name, 0
      name = f + "/" + fname + ".cpp"
      if os.path.isfile(name):
        return name, 0
      # if have a prefix, remove and try
      if prefix and fname.startswith(prefix):
        fname = fname[len(prefix):]
        name = f + "/" + fname + ".cc"
        if os.path.isfile(name):
          return name, 0
        name = f + "/" + fname + ".cpp"
        if os.path.isfile(name):
          return name, 0

  # if here, we still dont have a file match
  # if allowed to scan files, do that
  if scanfiles:
    #sys.stderr.write("Warning: Scaning for {}\n".format(fname))
    for f in paths:
      files = list(f + "/" + a for a in os.listdir(f))
      cc_files = fnmatch.filter(files, "*.cc")
      cpp_files = fnmatch.filter(files, "*.cpp")
      cxx_files = fnmatch.filter(files, "*.cxx")

      for fn in cc_files + cpp_files + cxx_files:
        line = find_defun_line_in_file(fn, fname)
        if line >= 0:
          #sys.stderr.write("Warning: Found function for {} in {} at {}\n".format(fname, fn, line))
          return fn, line
  
  return None, -1

def display_standalone_header():
  # make a file that doesnt need to be included in a texinfo file to
  # be valid
  print("@c mkfuncdocs output for a standalone function list")
  print("@include macros.texi")
  print("@ifnottex")
  print("@node Top")
  print("@top Function Documentation")
  print("Function documentation extracted from texinfo source in octave source files.")
  print("@contents")
  print("@end ifnottex")
  print("@node Function Reference")
  print("@chapter Function Reference")
  print("@cindex Function Reference")

def display_standalone_footer():
  print("@bye")

def display_func(name, ref, help):
  print ("@c -----------------------------------------")
  print ("@subsection {}".format(name))
  print ("@cindex {}".format(ref))
  for l in help:
    print ("{}".format(l))

def process (args):
  options = { 
    "verbose": False,
    "srcdir": [],
    "funcprefix": "",
    "ignore": [],
    "standalone": False,
    "allowscan": False
  }
  indexfile = ""

  for a in args:
    #print ("{}".format(a))
    c=a.split("=")
    key=c[0]

    if len(c) > 1:
      val=c[1]
    else:
      val=""

    if key == "--verbose":
      options["verbose"] = True;
    if key == "--standalone":
      options["standalone"] = True;
    elif key == "--allowscan":
      options["allowscan"] = True;
    elif key == "--src-dir":
      if val:
        options["srcdir"].append(val);
    elif key == "--ignore":
      if val:
        options["ignore"].append(val);
    elif key == "--func-prefix":
      if val:
        options["funcprefix"] = val;
    elif val == "":
      if indexfile == "":
        indexfile = key

  if indexfile == "":
    raise Exception("No index filename")

  if len(options["srcdir"]) == 0:
    options["srcdir"].append("inst")

  #print "options=", options
  if options['standalone']:
      display_standalone_header()

  idx = read_index(indexfile,  options["ignore"])
  for g in idx.groups:
    #print ("************ {}".format(g.name))
    g_name = texify_line(g.name)
    print ("@c ---------------------------------------------------")
    print ("@node {}".format(g_name))
    print ("@section {}".format(g_name))
    print ("@cindex {}".format(g_name))

    for f in sorted(g.functions):
      print ("@c {} {}".format(g_name, f))
      h = ""
      filename = ""
      path = ""
      if "@" in f:
        #print ("class func")
        path = f
        name = "@" + f
        ref = f.split("/")[-1]
        filename, lineno = find_func_file(path, options["srcdir"], options["funcprefix"])
      elif "." in f:
        path = f
        ref = f.split(".")[-1]
        name = f.split(".")[-1]
        filename, lineno = find_class_file(path, options["srcdir"])

        if not filename:
          parts = f.split('.')
          cnt  = 0
          path = ""
          for p in parts:
            if cnt < len(parts)-1:
              path = path + "/+"
            else:
              path = path + "/"
            path = path + p
            cnt = cnt + 1
          name = f;
          ref = parts[-1]
          filename, lineno = find_func_file(path, options["srcdir"], options["funcprefix"])

      elif "/" in f:
        path = f
        name = f
        ref = f.split("/")[-1]
        filename, lineno = find_func_file(path, options["srcdir"], options["funcprefix"])
      else:
        path = f
        name = f
        ref = f
        filename, lineno = find_func_file(path, options["srcdir"], options["funcprefix"], options['allowscan'])

      if not filename:
        sys.stderr.write("Warning: Cant find source file for {}\n".format(f))
      else:
        h = read_help (filename, lineno)

      if h:
        display_func (name, ref, h)

  if options['standalone']:
      display_standalone_footer()


def show_usage():
  print (sys.argv[0], "[options] indexfile")

if __name__ == "__main__":
  if len(sys.argv) > 1:
    status = process(sys.argv[1:])
    sys.exit(status)
  else:
    show_usage()
