function attenuation_factors = attenuationCT_to_511(Nx, Ny, Nz, KVP, varargin)
%attenuationCT_to_511 Scales the CT(HU)-values to corresponding 511 keV
%attenuation coefficients. Uses trilinear interpolation to scale the
%HU-values to 511 keV attenuation coefficients. Supports any X-ray tube
%potentials, but values outside of 80, 100, 120 and 140 kVp are
%interpolated. The calculated attenuation coefficients can also be saved on
%disk.
%
% Examples:
%   attenuation_factors = attenuationCT_to_511(Nx, Ny, Nz, KVP)
%   attenuation_factors = attenuationCT_to_511(Nx, Ny, Nz, KVP, attenuation_data, name)
% INPUTS:
%   Nx = Image dimension in x-axis
%   Ny = Image dimension in y-axis
%   Nz = Image dimension in z-axis
%   KVP = The CT kVp value (e.g. 100)
%   attenuation_data = The data that needs to be converted to 511 keV
%   attenuation values. If omitted or left empty, the user will be prompted
%   to load the data (MAT-file).
%   name = Examination name (optional). If specified, will be used to name
%   the MAT-file containing the attenuation images. If omitted or left
%   empty, no file will be saved.
%
% OUTPUTS:
%   attenuation_factors = 511 keV attenuation coefficients
% 
% See also attenuation122_to_511, create_atten_matrix_CT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ville-Veikko Wettenhovi, Samuli Summala
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Scale the attenuation coefficients %%%
%%% trilinear interpolation            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optargs = {[], []};
optargs(1:length(varargin)) = varargin;
reresized = optargs{1};
name = optargs{2};
% Values from: Accuracy of CT-based attenuation correction in PET/CT bone imaging
if KVP == 100
    a = [9.3 4 0.5]*10^(-5);
    b = [0.093 0.093 0.128];
elseif KVP == 80
    a = [9.3 3.28 0.41]*10^(-5);
    b = [0.093 0.093 0.122];
elseif KVP == 120
    a = [9.3 4.71 0.589]*10^(-5);
    b = [0.093 0.093 0.134];
elseif KVP == 140
    a = [9.3 5.59 0.698]*10^(-5);
    b = [0.093 0.093 0.142];
else
    warning('Unsupported kVp, interpolating initial values')
    a1 = [9.3 3.28 0.41]*10^(-5);
    b1 = [0.093 0.093 0.122];
    a2 = [9.3 4 0.5]*10^(-5);
    b2 = [0.093 0.093 0.128];
    a3 = [9.3 4.71 0.589]*10^(-5);
    b3 = [0.093 0.093 0.134];
    a4 = [9.3 5.59 0.698]*10^(-5);
    b4 = [0.093 0.093 0.142];
    aa = [a1; a2; a3; a4];
    bb = [b1; b2; b3; b4];
    c = [80;100;120;140];
    a = zeros(3,1);
    b = zeros(3,1);
    for kk = 1 : 3
        a(kk) = interp1(c,aa(:,kk),KVP,'spline');
        b(kk) = interp1(c,bb(:,kk),KVP,'spline');
    end
end
if isempty(reresized)
    [file, fpath] = uigetfile('*.mat','Select the CT images');
    if isequal(file, 0)
        error('No file was selected')
    end
    data = load([fpath file]);
    variables = fieldnames(data);
    reresized = double(data.(variables{1}));
    clear data
end
if size(reresized,3) ~= Nz && size(reresized,2) ~= Ny && size(reresized,1) ~= Nx
    error('The size of the input image does not match the input dimensions')
end
if size(reresized,3) ~= Nz
    reresized = reshape(reresized, Nx, Ny, Nz);
end
%risteyskohdat:
x(1,:) = [-1000,b(1)-1000*a(1)];
x(2,:) = [0,b(2)];
x(3,:) = [1000,b(2)+a(2)*1000];
x(4,:) = [3000,b(2)+a(2)*1000+a(3)*2000];
tarkkuus = 0.1;
inter = interp1(x(:,1),x(:,2),-1000:tarkkuus:3000,'linear');
vali = -1000:tarkkuus:3000;
% plot(vali,inter)
%asetetaan sopivasti
attenuation_factors = interp1(vali',inter',reresized(:),'linear','extrap');
attenuation_factors = reshape(attenuation_factors,Nx,Ny,Nz);  
if ~isempty(name)
    save(['511keV_attenuation_coefficients_for_' name '.mat'], 'attenuation_factors')
end