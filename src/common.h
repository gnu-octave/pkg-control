/*

Copyright (C) 2009-2016   Lukas F. Reichlin

This file is part of LTI Syncope.

LTI Syncope is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LTI Syncope is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LTI Syncope.  If not, see <http://www.gnu.org/licenses/>.

Common code for oct-files.

Author: Lukas Reichlin <lukas.reichlin@gmail.com>
Created: February 2012
Version: 0.2

*/

#ifndef COMMON_H
#define COMMON_H

#include <octave/f77-fcn.h>

#if defined (OCTAVE_HAVE_F77_INT_TYPE)
#  define TO_F77_INT(x) octave::to_f77_int (x)
#else
typedef octave_idx_type F77_INT;
#  define TO_F77_INT(x) (x)
#endif

F77_INT max (F77_INT a, F77_INT b);
F77_INT max (F77_INT a, F77_INT b, F77_INT c);
F77_INT max (F77_INT a, F77_INT b, F77_INT c, F77_INT d);
F77_INT max (F77_INT a, F77_INT b, F77_INT c, F77_INT d, F77_INT e);
F77_INT min (F77_INT a, F77_INT b);

void error_msg (const char name[], octave_idx_type index, octave_idx_type max, const char* msg[]);
void warning_msg (const char name[], octave_idx_type index, octave_idx_type max, const char* msg[]);
void warning_msg (const char name[], octave_idx_type index, octave_idx_type max, const char* msg[], octave_idx_type offset);

#endif
