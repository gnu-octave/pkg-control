#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([control], [3.2.0])
AC_CONFIG_SRCDIR([sl_ab01od.cc])
AC_CONFIG_HEADERS([config.h])

# Avoid warnings for redefining AH-generated preprocessor symbols of
# Octave.
AH_TOP([#include "undef-ah-octave.h"])

AC_CONFIG_MACRO_DIRS([m4])

# Checks for programs.
AC_CHECK_PROG(MKOCTFILE, mkoctfile, mkoctfile)
if test -z "$MKOCTFILE"; then
AC_MSG_ERROR([mkoctfile not found], 1);
fi
AC_CHECK_PROG(OCTAVE_CONFIG, octave-config, octave-config)
if test -z "$OCTAVE_CONFIG"; then
AC_MSG_ERROR([octave-config not found], 1);
fi

AC_PROG_CXX

# Start of checks for Octave features, preparations for checks.
OCTLIBDIR=${OCTLIBDIR:-`$OCTAVE_CONFIG -p OCTLIBDIR`}
## We need Octaves include path both with and without '/octave'
## appended. The path without '/octave' is needed to selectively test
## for Octave headers, like octave/....h. The path with '/octave' is
## needed since some Octave headers contain include directives for
## other Octave headers with <> instead of "".
OCTINCLUDEDIR=${OCTINCLUDEDIR:-`$MKOCTFILE -p INCFLAGS`}
AC_LANG_PUSH([C++])
TCXXFLAGS=$CXXFLAGS
TLDFLAGS=$LDFLAGS
TLIBS=$LIBS
TCPPFLAGS=$CPPFLAGS
LDFLAGS="-L$OCTLIBDIR $LDFLAGS"
LIBS="-loctinterp $LIBS"
# CXXFLAGS=
CPPFLAGS="$OCTINCLUDEDIR $CPPFLAGS"

## Simple symbol alternatives of different Octave versions.
OF_OCTAVE_LIST_ALT_SYMS([
[dnl
  [is_real_type],
  [isreal],
  [[octave_value ().isreal ();]],
  [OV_ISREAL],
  [],
  []
],

[dnl
  [is_cell],
  [iscell],
  [[octave_value ().iscell ();]],
  [OV_ISCELL],
  [],
  []
],

[dnl
  [is_object],
  [isobject],
  [[octave_value ().isobject ();]],
  [OV_ISOBJECT],
  [],
  []
],

[dnl
  [is_complex_type],
  [iscomplex],
  [[octave_value ().iscomplex ();]],
  [OV_ISCOMPLEX],
  [],
  []
],

[dnl
  [is_numeric_type],
  [isnumeric],
  [[octave_value ().isnumeric ();]],
  [OV_ISNUMERIC],
  [],
  []
]

],
   [oct-alt-includes.h])


LIBS=$TLIBS
LDFLAGS=$TLDFLAGS
CXXFLAGS=$TCXXFLAGS
CPPFLAGS=$TCPPFLAGS
AC_LANG_POP([C++])
# End of checks for Octave features.

AC_OUTPUT
