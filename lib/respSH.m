function [h1, h2, h3] = respSH(Depth, Vs, Rho, Fox, f, th)
% RESPSH Returns seismic response of damped soil layers to incident SH-wave
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Returns seismic response of damped soil layers on rock to incident SH-wave
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Miroslav HALLO
% ETH Zürich, Swiss Seismological Service
% E-mail: miroslav.hallo@sed.ethz.ch
% Revision 2019/06: The first version of the function
% Tested in Matlab R2018b, R2019b, R2021a, 2025b
% Method by Kramer, S.L. (1996): Geotechnical Earthquake Engineering, 
%           Prentice Hall, pp. 268-270.
%
% Copyright (C) 2019 Swiss Seismological Service, ETH Zurich
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
% Fox - Layer damping ratio
% f - frequency discretitation [Hz]
% th - Wave incidence angle [deg] (0-60, beware higher angles - not realistic)
%
% OUTPUT:
% h1 - "surface to incident SH" transfer function
% h2 - "surface to borehole" transfer function
% h3 - "surface to rock outcrop" transfer function
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare variables
nL = length(Depth); % Number of layers
nf = length(f); % Number of discrete freqencies

Thick = zeros(1,nL); % Thickness of layers
Thick(1:end-1) = Depth(2:end) - Depth(1:end-1);
Thick(end) = Inf;

VsS = Vs.*(1+1i*Fox); % Complex shear wave velocity

Hslow = sind(th)/Vs(end); % Horizontal inc. shear wave slowness
Vslow = sqrt(( 1./(VsS.^2)) - (Hslow^2) ); % Vertical inc. shear wave slowness

% Complex impedance ratio
VsSimp = (VsS.^2) .* Vslow;
alphaS = (Rho(1:end-1).*VsSimp(1:end-1)) ./ (Rho(2:end).*VsSimp(2:end)); 

A = complex(ones(nL,nf));
B = complex(ones(nL,nf));

% Recursive loop over layers
for m = 1:nL-1
    % Complex wave number x thickness
    ksH = 2*pi*f*Thick(m)*Vslow(m);
    % Amplitudes of up-going and down-going waves
    A(m+1,1:nf) = 0.5*A(m,1:nf)*(1+alphaS(m)).*exp(1i*ksH) + 0.5*B(m,1:nf)*(1-alphaS(m)).*exp(-1i*ksH);
    B(m+1,1:nf) = 0.5*A(m,1:nf)*(1-alphaS(m)).*exp(1i*ksH) + 0.5*B(m,1:nf)*(1+alphaS(m)).*exp(-1i*ksH);
end

% Transfer function "surface to incident SH"  (A at nL(top) -> surface)
h1(1:nf) = 2 ./ A(nL,:);

% Transfer function "surface to borehole"     (A+B at nL(top) -> surface)
h2(1:nf) = ( A(1,:) + B(1,:) ) ./ (  A(nL,:) + B(nL,:) );

% Transfer function "surface to rock outcrop" (A at nL(surface) -> surface)
h3(1:nf) = 1 ./ A(nL,:);

end

