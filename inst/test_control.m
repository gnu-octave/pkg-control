## Copyright (C) 2009-2016   Lukas F. Reichlin
##
## This file is part of LTI Syncope.
##
## LTI Syncope is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## LTI Syncope is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with LTI Syncope.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Script File} {} test_control
## Execute all available tests at once.
##
## The Octave control package uses the 
## Uses @url{https://github.com/SLICOT/SLICOT-Reference, SLICOT library},
## Copyright (c) 1996-2025, SLICOT, available under the BSD 3-Clause
## (@url{https://github.com/SLICOT/SLICOT-Reference/blob/main/LICENSE,  License and Disclaimer}).
## @acronym{SLICOT} needs @acronym{BLAS} and @acronym{LAPACK} libraries which are also prerequisites
## for Octave itself.
## In case of failing tests, it is highly recommended to use
## @url{http://www.netlib.org/blas/, Netlib's reference @acronym{BLAS}} and
## @url{http://www.netlib.org/lapack/, @acronym{LAPACK}}
## for building Octave.  Using @acronym{ATLAS} may lead to sign changes
## in some entries of the state-space matrices.
## In general, these sign changes are not 'wrong' and can be regarded as
## the result of state transformations.  Such state transformations
## (but not input/output transformations) have no influence on the
## input-output behaviour of the system.  For better numerics,
## the control package uses such transformations by default when
## calculating the frequency responses and a few other things.
## However, arguments like the Hankel singular Values (@acronym{HSV}) must not change.
##
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: May 2010
## Version: 0.6
##
## Modified 01.01.2024 by Torsten Lilge <ttl-octave@mailbox.org>
## Use pkg test control instead of calling all tests separatly

pkg test control
