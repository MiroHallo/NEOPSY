function [amin,amax,fticks] = fscale(fmin,fmax)
% FSCALE Return frequency axis limits and intelligent ticks
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Miroslav HALLO
% ETH Zürich, Swiss Seismological Service
% E-mail: miroslav.hallo@sed.ethz.ch
% Revision 4/2021: The first version
%
% Copyright (C) 2021  Swiss Seismological Service, ETH Zurich
%
% This program is published under the GNU General Public License (GNU GPL).
%
% This program is free software: you can modify it and/or redistribute it
% or any derivative version under the terms of the GNU General Public
% License as published by the Free Software Foundation, either version 3
% of the License, or (at your option) any later version.
%
% This code is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
% and don't remove their names from the code.
%
% You should have received copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
% fmin - Min frequency [Hz]
% fmax - Max frequency [Hz]
%
% OUTPUT:
% amin - Min anchored frequency [Hz]
% amax - Max anchored frequency [Hz]
% fticks - Intelligent ticks [Hz]
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants

% Define stable anchor points for axis limits
stable_limits = [0.05 0.1, 0.2, 0.5, 1, 5, 10, 20, 30, 50, 100];

% Find closest lower and upper anchor points
amin = stable_limits(find(stable_limits <= fmin, 1, 'last'));
amax = stable_limits(find(stable_limits >= fmax, 1, 'first'));

% If outside the range
if isempty(amin), amin = min(stable_limits); end
if isempty(amax), amax = max(stable_limits); end

% Generate all ticks
decades = 10.^(-2:1);
major_prio = [1, 2, 5]';
minor_prio = [3, 4, 6, 8]';
detail_prio = [1.5, 2.5, 7, 9]';

ticks_major  = sort(reshape(major_prio  * decades, 1, []));
ticks_minor  = sort(reshape(minor_prio  * decades, 1, []));
ticks_detail = sort(reshape(detail_prio * decades, 1, []));

% Select tick set
ratio = amax/amin;
if ratio > 10
    all_ticks = ticks_major;
elseif ratio > 3
    all_ticks = [ticks_major, ticks_minor];
else
    all_ticks = [ticks_major, ticks_minor, ticks_detail];
end
all_ticks = sort(all_ticks);

% Select only ticks that fall within limits
fticks = all_ticks(all_ticks >= amin & all_ticks <= amax);

end

