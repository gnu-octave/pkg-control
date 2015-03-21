/*

Copyright (C) 2009-2015   Lukas F. Reichlin

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

octave_idx_type max (octave_idx_type a, octave_idx_type b);
octave_idx_type max (octave_idx_type a, octave_idx_type b, octave_idx_type c);
octave_idx_type max (octave_idx_type a, octave_idx_type b, octave_idx_type c, octave_idx_type d);
octave_idx_type max (octave_idx_type a, octave_idx_type b, octave_idx_type c, octave_idx_type d, octave_idx_type e);
octave_idx_type min (octave_idx_type a, octave_idx_type b);
void error_msg (const char name[], octave_idx_type index, octave_idx_type max, const char* msg[]);
void warning_msg (const char name[], octave_idx_type index, octave_idx_type max, const char* msg[]);
void warning_msg (const char name[], octave_idx_type index, octave_idx_type max, const char* msg[], octave_idx_type offset);

#endif
