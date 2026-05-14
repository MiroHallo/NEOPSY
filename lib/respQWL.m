function [z_f,Vs_f,Rho_f,A,z2_f,Vs2_f,IC_f] = respQWL(Depth,Vs,Rho,f)
% RESPQWL Returns amplification by quarter-wavelength average velocity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computing the amplification by quarter-wavelength average velocity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Miroslav HALLO
% ETH Zürich, Swiss Seismological Service
% E-mail: miroslav.hallo@sed.ethz.ch
% Revision 6/2019: The first version of the function.
% Method by:
% Boore, D.M. (2003): Simulation of Ground Motion Using the Stochastic Method,
%       Pure and Applied Geophysics, 160(3-4), 635-676.
% Poggi, V., Edwards, B., Fah, D. (2011): Derivation of a Reference Shear-Wave
%       Velocity Model from Empirical Site Amplification, Bulletin of the 
%       Seismological Society of America, 101(1), 258-274.
% Poggi, V., Edwards, B., Fah, D. (2012): Characterizing the Vertical-to-Horizontal
%       Ratio of Ground Motion at Soft-Sediment Sites, Bulletin of the 
%       Seismological Society of America, 102(6), 2741-2756.
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
% Rho - Density [kg/m3]
% f - frequency discretitation [Hz]
%
% OUTPUT:
% z_f - QWL depth
% Vs_f - QWL velocity
% Rho_f - QWL density
% A - The amplification factor
% z2_f - 2nd QWL depth
% Vs2_f - 2nd QWL velocity
% IC_f - QWL Impedance contrast
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VsC = Vs(end); % Reference (source) shear wave velocity [m/s]
RhoC = Rho(end); % Reference (source) density [kg/m3]
z = logspace(-2,5,10000); % Depth search discretitation [m]

% -----------------------------------------------------------------
% Prepare variables
nL = length(Depth);   % number of layers
nz = length(z);   % number of discrete depths
nf = length(f);   % number of discrete freqencies
Thick = zeros(1,nL); % thickness of layers
Thick(1:end-1) = Depth(2:end) - Depth(1:end-1);
Thick(end) = Inf;

% Find layer indexes
IZ = zeros(1,nz)+nL;
IZ_tmp = 1;
ni = 0;
while 1
    ni = ni+1;
    if (IZ_tmp == nL) || (ni>nz)
        break
    elseif Depth(IZ_tmp+1)>=z(ni)
        IZ(ni) = IZ_tmp;
    else
        IZ_tmp = IZ_tmp + 1;
        ni = ni-1;
    end
end

% Compute vertical travel-times
T0 = [0 cumsum(Thick./Vs)];
T = T0(IZ) + abs(z - Depth(IZ)) ./ Vs(IZ);

% Find quarter-wavelength depth and velocity
VsQWL = z./T;
minI = zeros(1,nf);

for i = 1:nf
    [~,minI(i)] = min( abs(z - VsQWL/(4*f(i))) );
end

z_f = z(minI);
Vs_f = VsQWL(minI);
IZ_f = IZ(minI);

% Find quarter-wavelength density
Rho0 = [0 cumsum(Rho.*Thick)];
Rho_f = (Rho0(IZ_f) + abs(z_f - Depth(IZ_f)) .* Rho(IZ_f)) ./ z_f;

% Amplification
A = sqrt((VsC*RhoC)./(Vs_f.*Rho_f)); 

% Find quarter-wavelength Impedance contrast
z2_f = zeros(1,nf);
Vs2_f = zeros(1,nf);
for i = 1:nf
    T_tmp = T(minI(i):end) - T(minI(i));
    z_tmp = z(minI(i):end) - z(minI(i));
    VsQWL_tmp = z_tmp./T_tmp;
    
    [~,minC] = min(abs(z_tmp - VsQWL_tmp/(4*f(i))));
    
    z2_f(i) = z(minI(i)-1+minC);
    Vs2_f(i) = VsQWL_tmp(minC);
end
IC_f = Vs_f./Vs2_f;

end

