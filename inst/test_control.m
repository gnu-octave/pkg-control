## Copyright (C) 2010, 2011   Lukas F. Reichlin
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
## @deftypefn {Script File} test_control
## Execute all available tests at once.
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: May 2010
## Version: 0.2.1

## test collection
test ltimodels

## robust control
test hinfsyn
test h2syn
test ncfsyn

## ARE solvers
test care
test dare
test kalman

## Lyapunov
test lyap
test dlyap
test gram
test covar
test lyapchol
## test dlyapchol  # TODO: add tests

## various oct-files
test place
test hsvd

## various m-files
test margin
test sigma
test initial
test ctrb
test obsv
test issample