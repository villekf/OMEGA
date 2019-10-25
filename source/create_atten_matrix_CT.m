function attenuation_factors = create_atten_matrix_CT(Nx, Ny, Nz, axial_fov, FOVa_x, varargin)
%CREATE_ATTEN_MATRIX_CT Loads DICOM CT images (.ima-files) and scales the
% CT(HU)-values to corresponding 511 keV attenuation coefficients. Uses
% trilinear interpolation to scale the HU-values to 511 keV attenuation
% coefficients. Supports X-ray tube potentials of 80, 100, 120 and 140 kVp.
% The user will be prompted for the location (folder) of the CT images. It
% is assumed that the PET FOV will be centered in the center of the (larger)
% CT FOV. The calculated attenuation coefficients are also saved on disk if
% a name for the examination is provided.
%
% This code requires that the image width information in the DICOM image
% is named as Width, image height as Height, pixel spacing as PixelSpacing
% and slice thickness as SliceThickness. 
% EXPERIMENTAL CODE.
%
% Example:
%   attenuation_factors = attenuationCT_to_511(options)
% INPUTS:
%   Nx = Image dimension in x-axis
%   Ny = y-axis
%   Nz = z-axis
%   axial_fov = FOV size (side length) in axial direction
%   FOV = FOV size (side length) in transaxial direction
%   KVP = The CT kVp value (e.g. 100). This parameter is optional. If the
%   DICOM image contains info named as KVP this can be left empty.
%   name = Examination name (optional). If specified, will be used to name
%   the MAT-file containing the attenuation images. If omitted or left
%   empty, no file will be saved.
%
% OUTPUTS:
%   attenuation_factors = 511 keV attenuation coefficients
% 
% See also attenuation122_to_511, attenuationCT_to_511


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi, Samuli Summala
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
optargs = {[], []};
optargs(1:length(varargin)) = varargin;
KVP = optargs{1};
name = optargs{2};
fpath = uigetdir([], 'Select folder with .ima CT-images');
if isequal(fpath, 0)
    error('No folder was selected')
end
fnames = dir([fpath '/*.IMA']);
numfids = length(fnames);
info = dicominfo([fnames(1).folder '/' fnames(1).name]);
W = info.Width;    %width of image
H = info.Height;    %height of image
if isempty(KVP)
    KVP = info.KVP;
end
atten_image = zeros(H,W,length(numfids));
for k = 1:numfids
  atten_image(:,:,k) = dicomread([fnames(k).folder '/' fnames(k).name]);
end
pix_space = info.PixelSpacing;
slice_thick = info.SliceThickness;
Pet_pixels = Nx;
Nslices = Nz;
slice_lkm_ct = floor(axial_fov/slice_thick) + 2*(ceil(((numfids*slice_thick - axial_fov)/2)/slice_thick) - floor(((numfids*slice_thick - axial_fov)/2)/slice_thick));
atten_image = atten_image(:,:,slice_lkm_ct+1:222-slice_lkm_ct)-1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale the CT images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resized = zeros(W,H,Nslices);
for i=1:W
    resized(:,i,:) = reshape(imresize(reshape(atten_image(:,i,:),W,slice_lkm_ct),[W Nslices]),W,1,Nslices);
end
pixel_lkm_ct = ceil(double(W) - ((double(W)*pix_space(1)) - FOVa_x) / pix_space(1));
dd = floor((double(W)-pixel_lkm_ct)/2);
resized = resized(dd+2:end-dd,dd+2:end-dd,:);
reresized = zeros(Pet_pixels,Pet_pixels,Nslices);
for i=1:Nslices
    reresized(:,:,i) = reshape(imresize(reshape(resized(:,:,i),pixel_lkm_ct,pixel_lkm_ct),[Pet_pixels Pet_pixels]),Pet_pixels,Pet_pixels,1);
end
%%

attenuation_factors = attenuationCT_to_511(Nx, Ny, Nz, KVP, reresized, name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Scale the attenuation coefficients %%%
%%% trilinear interpolation            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Values from: Accuracy of CT-based attenuation correction in PET/CT bone imaging
% if KVP == 100
%     a = [9.3 4 0.5]*10^(-5);
%     b = [0.093 0.093 0.128];
% elseif KVP == 80
%     a = [9.3 3.28 0.41]*10^(-5);
%     b = [0.093 0.093 0.122];
% elseif KVP == 120
%     a = [9.3 4.71 0.589]*10^(-5);
%     b = [0.093 0.093 0.134];
% elseif KVP == 140
%     a = [9.3 5.59 0.698]*10^(-5);
%     b = [0.093 0.093 0.142];
% else
%     error('Unsupported kVp')
% end
% %risteyskohdat:
% x(1,:) = [-1000,b(1)-1000*a(1)];
% x(2,:) = [0,b(2)];
% x(3,:) = [1000,b(2)+a(2)*1000];
% x(4,:) = [3000,b(2)+a(2)*1000+a(3)*2000];
% tarkkuus = 0.1;
% inter = interp1(x(:,1),x(:,2),-1000:tarkkuus:3000,'linear');
% vali = -1000:tarkkuus:3000;
% % plot(vali,inter)
% %asetetaan sopivasti
% attenuation_factors = zeros(size(reresized));
% for ii=1:Nslices
% %     ii
%     temppi=reresized(:,:,ii);
%     temppi=temppi(:);
%     tic
%     for ll=1:Pet_pixels^2
%         apu=temppi(ll);
%         [~,idx]=min(abs(vali-apu));
%         apu=inter(idx);
%         temppi(ll)=apu;
%     end
%     toc
%     attenuation_factors(:,:,ii)=reshape(temppi,Pet_pixels,Pet_pixels);  
% end
% % imagesc(attenuation_factors(:,:,55))
% save(['attenuation_coefficients_for_' name '.mat'], attenuation_factors)