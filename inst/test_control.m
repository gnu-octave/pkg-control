## Copyright (C) 2010   Lukas F. Reichlin
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
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Script File} test_control
## Execute all available tests at once.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: May 2010
## Version: 0.1

test ltimodels
test hinfsyn
test h2syn
test place
test care

test kalman
test lyap
test dlyap
test gram
test covar

test margin
test sigma
test initial
test hsvd