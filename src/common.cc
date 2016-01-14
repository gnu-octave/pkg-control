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
Created: April 2010
Version: 0.4

*/


#include <octave/oct.h>

octave_idx_type max (octave_idx_type a, octave_idx_type b)
{
    if (a > b)
        return a;
    else
        return b;
}

octave_idx_type max (octave_idx_type a, octave_idx_type b, octave_idx_type c)
{
    return max (max (a, b), c);
}

octave_idx_type max (octave_idx_type a, octave_idx_type b, octave_idx_type c, octave_idx_type d)
{    
    return max (max (a, b), max (c, d));
}

octave_idx_type max (octave_idx_type a, octave_idx_type b, octave_idx_type c, octave_idx_type d, octave_idx_type e)
{
    return max (max (a, b, c, d), e);
}

octave_idx_type min (octave_idx_type a, octave_idx_type b)
{
    if (a < b)
        return a;
    else
        return b;
}

void error_msg (const char name[], octave_idx_type index, octave_idx_type max, const char* msg[])
{
    if (index == 0)
        return;

    if (index < 0)
        error ("%s: the %d-th argument had an invalid value", name, index);
    else if (index <= max)
        error ("%s: %s", name, msg[index]);
    else
        error ("%s: unknown error, info = %d", name, index);
}

void warning_msg (const char name[], octave_idx_type index, octave_idx_type max, const char* msg[])
{
    if (index == 0)
        return;

    if (index > 0 && index <= max)
        warning ("%s: %s", name, msg[index]);
    else
        warning ("%s: unknown warning, iwarn = %d", name, index);
}

void warning_msg (const char name[], octave_idx_type index, octave_idx_type max, const char* msg[], octave_idx_type offset)
{
    if (index == 0)
        return;

    if (index > 0 && index <= max)
        warning ("%s: %s", name, msg[index]);
    else if (index > offset)
        warning ("%s: %d+%d: %d %s", name, offset, index-offset, index-offset, msg[max+1]);
    else
        warning ("%s: unknown warning, iwarn = %d", name, index);
}
