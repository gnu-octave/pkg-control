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

#include "common.h"

F77_INT max (F77_INT a, F77_INT b)
{
    if (a > b)
        return a;
    else
        return b;
}

F77_INT max (F77_INT a, F77_INT b, F77_INT c)
{
    return max (max (a, b), c);
}

F77_INT max (F77_INT a, F77_INT b, F77_INT c, F77_INT d)
{    
    return max (max (a, b), max (c, d));
}

F77_INT max (F77_INT a, F77_INT b, F77_INT c, F77_INT d, F77_INT e)
{
    return max (max (a, b, c, d), e);
}

F77_INT min (F77_INT a, F77_INT b)
{
    if (a < b)
        return a;
    else
        return b;
}

void error_msg (const char name[], int index, octave_idx_type max, const char* msg[])
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

void warning_msg (const char name[], int index, octave_idx_type max, const char* msg[])
{
    if (index == 0)
        return;

    if (index > 0 && index <= max)
        warning ("%s: %s", name, msg[index]);
    else
        warning ("%s: unknown warning, iwarn = %d", name, index);
}

void warning_msg (const char name[], int index, octave_idx_type max, const char* msg[], int offset)
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
