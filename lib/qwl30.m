function [f30,Vs30] = qwl30(Depth,Vs,d30)
% QWL30 Computing the Vs30 of layered velocity model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Miroslav HALLO
% ETH Zürich, Swiss Seismological Service
% E-mail: miroslav.hallo@sed.ethz.ch
% Revision 10/2019: The first version of the function.
%
% Copyright (C) 2019  Swiss Seismological Service, ETH Zurich
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
% Depth -  Top depth of layers [m]
% Vs - Shear wave velocity [m/s]
% d30 - depth of interest (30) [m]
%
% OUTPUT:
% f30 - Frequency of QWL representation of the VS30
% Vs30 - VS30 velocity [m/s]
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare variables
nL = length(Depth);   % number of layers
Thick = zeros(1,nL); % thickness of layers
Thick(1:end-1) = Depth(2:end) - Depth(1:end-1);
Thick(end) = Inf;

% Find layer indexes
IZ = nL;
IZ_tmp = 1;
ni = 0;
while 1
    ni = ni+1;
    if (IZ_tmp == nL) || (ni>1)
        break
    elseif Depth(IZ_tmp+1)>=d30
        IZ(ni) = IZ_tmp;
    else
        IZ_tmp = IZ_tmp + 1;
        ni = ni-1;
    end
end

% Compute vertical travel-times
T0 = [0 cumsum(Thick./Vs)];
T = T0(IZ) + abs(d30 - Depth(IZ)) ./ Vs(IZ);

% Find the vs30 frequency and velocity
f30 = 1/(4*T);
Vs30 = d30 / T;

end

