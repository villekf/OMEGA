function img = attenuation122_to_511(img)
%ATTENUATION122_TO_511 Scales the 122 keV attenuation coefficients to 511
%keV. This code assumes that the attenuation coefficients are 1/cm. First
%calculate the tabulated values for various tissues and elements for both
%122 and 511 keV cases. Scale the values such that the peak is at the soft
%tissue level (ignore air). Set small values to air. Interpolate the
%densities corresponding to the attenuation coefficients obtained for 122
%keV. Interpolate the corresponding attenuation coefficients for 511 keV
%based on the interpolated densities. The calculated attenuation
%coefficients are also saved on disk.
%
% Example:
%   img = attenuationCT_to_511(img)
% INPUTS:
%   img = Attenuation image containing the attenuation coefficients for 122
%   keV
%
% OUTPUTS:
%   img = Attenuation image containing the attenuation coefficients for 511
%   keV
% 
% See also attenuationCT_to_511, create_atten_matrix_CT
%
% References
%   https://www.nist.gov/pml/x-ray-mass-attenuation-coefficients


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi
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

[Nx,Ny,Nz] = size(img);

bone122 = (1.851e-1+1.48e-1)/2;
bone122 = bone122 * 1.92;

bone511 = (9.022e-2 + (8.332e-2 - 9.022e-2)*0.11);
bone511 = bone511 * 1.92;

tissue122 = (1.693e-1 + 1.492e-1)/2;
tissue122 = tissue122 * 1.06;

tissue511 = (9.598e-2 + (8.873e-2 - 9.598e-2)*0.11);
tissue511 = tissue511 * 1.06;

adtissue122 = (1.688e-1 + 1.5e-1)/2;
adtissue122 = adtissue122 * 9.5e-1;

adtissue511 = (9.696e-2 + (8.965e-2 - 9.696e-2)*0.11);
adtissue511 = adtissue511 * 9.5e-1;

water122 = (1.707e-1 + 1.505e-1)/2;

water511 = (9.687e-2 + (8.956e-2 - 9.687e-2)*0.11);

muscle122 = (1.693e-1 + 1.492e-1)/2;
muscle122 = muscle122 * 1.05;

muscle511 = (9.598e-2 + (8.874e-2 - 9.598e-2)*0.11);
muscle511 = muscle511 * 1.05;

plastic122 = (1.793e-1 + 1.482e-1)/2;
plastic122 = plastic122 * 1.45;

plastic511 = 9.227e-2 + (8.525e-2 - 9.227e-2)*0.11;
plastic511 = plastic511 * 1.45;

air122 = (1.541e-1 + 1.356e-1)/2;
air122 = air122 * 1.205e-3;

air511 = (8.712e-2 + (8.055e-2 - 8.712e-2)*0.11);
air511 = air511 * 1.205e-3;

alu122 = (1.704e-1 + 1.378e-1)/2;
alu122 = alu122 * 2.699;

alu511 = (8.445E-02 + (7.802E-02 - 8.445E-02)*0.11);
alu511 = alu511 * 2.699;

[apu1,apu2] = hist(img(img > adtissue122/2));
[~,locmax] = max(apu1);
kerroin = tissue122/apu2(locmax);

img = img * kerroin;
img(img < adtissue122/100) = air122;

X = ([air122, adtissue122,water122,muscle122,tissue122,plastic122,bone122,alu122]);

V = [1.205e-3, 9.5e-1,1, 1.05, 1.06, 1.45, 1.92, 2.699];

Vq = interp1(X,V,img(:),'linear');

V = ([air511, adtissue511,water511,muscle511,tissue511,plastic511,bone511,alu511]);

X = [1.205e-3, 9.5e-1,1, 1.05, 1.06, 1.45, 1.92, 2.699];

Vq = interp1(X,V,Vq,'linear');

Vq(Vq < 1.205e-3) = 1.205e-3;

img = reshape(Vq,Nx,Ny,Nz);