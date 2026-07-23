## Copyright (C) 2026        Prateek Ganguli
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

## Author: Prateek Ganguli <prateek.ganguli@gmail.com>
## Created: July 2026
## Version: 0.1

## Extract the internal-delay port matrices (B2, C2, D12, D21, D22) and the
## per-port delay tau (whole samples) from a discrete InternalDelay ss model,
## reusing @ss/__ss_ext_build__ (the same helper c2d/d2c use) for data
## extraction only.  Shared by __time_response__.m and lsim.m.
function [B2, C2, D12, D21, D22, tau] = __internaldelay_ports__ (sys)
  [ext, nu, ny] = __ss_ext_build__ (sys);
  [~, Be, Ce, De] = ssdata (ext);
  B2  = Be(:, nu+1:end);
  C2  = Ce(ny+1:end, :);
  D12 = De(1:ny, nu+1:end);
  D21 = De(ny+1:end, 1:nu);
  D22 = De(ny+1:end, nu+1:end);
  tau = get (sys, "internaldelay")(:);
endfunction
