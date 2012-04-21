## Copyright (C) 2010, 2011, 2012   Lukas F. Reichlin
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
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: May 2010
## Version: 0.4

## test collection
test ltimodels

## LTI methods
test @lti/c2d
test @lti/d2c
test @lti/feedback
test @lti/horzcat
test @lti/inv
test @lti/minreal
test @lti/mtimes
test @lti/norm
test @lti/plus
test @lti/prescale
test @lti/sminreal
test @lti/zero

## robust control
test h2syn
test hinfsyn
test ncfsyn

## ARE solvers
test care
test dare
test kalman

## Lyapunov
test covar
test dlyap
## test dlyapchol  # TODO: add tests
test gram
test lyap
test lyapchol

## model order reduction
test bstmodred
test btamodred
test hnamodred
## test spamodred  # TODO: create test case

## controller order reduction
test btaconred
test cfconred
test fwcfconred
## test spaconred  # TODO: create test case

## identification
test fitfrd

## various oct-files
test ctrbf
test hsvd
test place

## various m-files
test ctrb
test filt
test initial
test issample
test margin
test obsv
test sigma
