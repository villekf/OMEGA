function attenuation_factors = create_atten_matrix_CT(options, varargin)
%CREATE_ATTEN_MATRIX_CT Loads DICOM CT images (.ima-files) and scales the
% CT(HU)-values to corresponding 511 keV attenuation coefficients. Uses
% trilinear interpolation to scale the HU-values to 511 keV attenuation
% coefficients. Supports X-ray tube potentials of 80, 100, 120 and 140 kVp.
% The user will be prompted for the location (folder) of the CT images. It
% is assumed that the PET FOV will be centered in the center of the CT FOV.
% The CT images are resampled to the PET image grid (Nx x Ny, FOVa_x x
% FOVa_y) with exact area weighting. If the PET FOV is larger than the CT
% FOV, the part outside the CT is filled with air. The calculated
% attenuation coefficients are also saved on disk if a name for the
% examination is provided.
%
% If options.scatter_correction is true, the narrow-beam water attenuation
% coefficient 0.096 cm^-1 is used in the conversion (the scattered events
% are modelled separately). Otherwise the effective value 0.093 cm^-1 is
% used, which partly compensates for the scatter remaining in the data. See
% attenuationCT_to_511 for details.
%
% The stored pixel values are converted to HU slice by slice with the
% DICOM RescaleSlope and RescaleIntercept (HU = slope * value + intercept).
% If a file contains no rescale information, the fixed offset 1000 is
% subtracted instead (HU = value - 1000) and a warning is shown once. HU
% values below -1000 (air) are set to -1000.
%
% This code requires that the image width information in the DICOM image
% is named as Width, image height as Height, pixel spacing as PixelSpacing
% and slice thickness as SliceThickness. 
% EXPERIMENTAL CODE.
%
% Example:
%   attenuation_factors = create_atten_matrix_CT(options, KVP, name)
% INPUTS:
%   options = The options struct (Nx, Ny, Nz, axial_fov, FOVa_x, FOVa_y,
%   attenuation_datafile and, optionally, scatter_correction and verbose
%   are used)
%   Nx = Image dimension in x-axis
%   Ny = y-axis
%   Nz = z-axis
%   axial_fov = FOV size (side length) in axial direction
%   FOVa_x = FOV size (side length) in transaxial direction (x-axis, i.e.
%   the first dimension of the CT image)
%   FOVa_y = FOV size in y-axis (optional field in options, FOVa_x is used if
%   not present or empty)
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
% Copyright (C) 2019 Ville-Veikko Wettenhovi, Samuli Summala
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
Nx = options.Nx;
Ny = options.Ny;
Nz = options.Nz;
axial_fov = options.axial_fov;
FOVa_x = options.FOVa_x;
if isfield(options, 'FOVa_y') && ~isempty(options.FOVa_y)
    FOVa_y = options.FOVa_y;
else
    FOVa_y = FOVa_x;
end
optargs = {[], []};
optargs(1:length(varargin)) = varargin;
KVP = optargs{1};
name = optargs{2};
if isempty(options.attenuation_datafile) || exist(options.attenuation_datafile,'dir') ~= 7
    fpath = uigetdir([], 'Select folder with .ima CT-images');
    if isequal(fpath, 0)
        error('No folder was selected')
    end
else
    fpath = options.attenuation_datafile;
end
fnames = dir([fpath '/*.IMA']);
numfids = length(fnames);
info = dicominfo([fnames(1).folder '/' fnames(1).name]);
% info2 = dicominfo([fnames(2).folder '/' fnames(2).name]);
W = info.Width;    %width of image
H = info.Height;    %height of image
if isempty(KVP)
    KVP = info.KVP;
end
atten_image = zeros(H,W,numfids);
noRescale = false;
for k = 1:1:numfids
  [atten_image(:,:,k), hasRescale] = readCTSliceHU([fnames(k).folder '/' fnames(k).name]);
  if ~hasRescale && ~noRescale
      noRescale = true;
      warning('create_atten_matrix_CT:noRescale', ...
          ['DICOM RescaleSlope/RescaleIntercept not found, using the offset -1000 to convert the CT values to HU (file ' fnames(k).name ')'])
  end
end
pix_space = info.PixelSpacing;
slice_thick = info.SliceThickness;
Nslices = Nz;
slice_lkm_alku = floor(((numfids * slice_thick - axial_fov) / 2) / slice_thick);
% atten_image = atten_image(:,:,slice_lkm_ct+1:end-slice_lkm_ct)-1000;
atten_image = atten_image(:,:,1 + slice_lkm_alku: end - slice_lkm_alku);
slice_lkm_ct = size(atten_image,3);
% atten_image = atten_image(:,:,1:slice_lkm_ct)-1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale the CT images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Axial scaling to Nz slices
resized = zeros(H,W,Nslices);
for i=1:W
    resized(:,i,:) = reshape(imresize(reshape(atten_image(:,i,:),H,slice_lkm_ct),[H Nslices]),H,1,Nslices);
end
% Transaxial resampling to the PET image grid by exact area weighting. The
% first dimension of the CT image is mapped to Nx pixels covering FOVa_x,
% the second to Ny pixels covering FOVa_y. Both grids are centered at the
% image center. The parts of the PET FOV outside the CT image are filled
% with air (-1000 HU).
if isscalar(pix_space)
    pix_space = [pix_space; pix_space];
end
Mx = overlapMatrix(H, pix_space(1), Nx, FOVa_x / Nx);
My = overlapMatrix(W, pix_space(2), Ny, FOVa_y / Ny);
coverage = min(max(full(Mx * ones(H,1)) * full(My * ones(W,1))', 0), 1);
airPart = (1 - coverage) * (-1000);
reresized = zeros(Nx,Ny,Nslices);
for i=1:Nslices
    reresized(:,:,i) = Mx * resized(:,:,i) * My' + airPart;
end
% The axial interpolation may undershoot slightly below air
reresized = max(reresized, -1000);
% reresized = flip(reresized,3);
%%

% Narrow-beam water attenuation (0.096 cm^-1) with scatter correction,
% effective value (0.093 cm^-1) otherwise
if isfield(options, 'scatter_correction') && options.scatter_correction
    waterMu = 0.096;
    if isfield(options, 'verbose') && options.verbose > 0
        disp('Scatter correction selected: using narrow-beam water attenuation 0.096 cm^-1 for the CT conversion')
    end
else
    waterMu = 0.093;
end
attenuation_factors = attenuationCT_to_511(Nx, Ny, Nz, KVP, reresized, name, waterMu);
attenuation_factors = attenuation_factors / 10;
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

function [hu, hasRescale] = readCTSliceHU(fname)
%READCTSLICEHU Reads a DICOM CT slice and converts it to HU
% HU = RescaleSlope * value + RescaleIntercept. If either field is missing
% or empty, HU = value - 1000 is used and hasRescale is false.
inf0 = dicominfo(fname);
hasRescale = isfield(inf0, 'RescaleSlope') && isfield(inf0, 'RescaleIntercept') ...
    && ~isempty(inf0.RescaleSlope) && ~isempty(inf0.RescaleIntercept);
if hasRescale
    hu = double(dicomread(fname)) * double(inf0.RescaleSlope) + double(inf0.RescaleIntercept);
else
    hu = double(dicomread(fname)) - 1000;
end
% Values below air (e.g. -1024 HU outside the reconstruction circle in
% Siemens CT) are not physical and would give negative attenuation
% coefficients in attenuationCT_to_511, which extrapolates linearly
hu = max(hu, -1000);

function M = overlapMatrix(nIn, psIn, nOut, psOut)
%OVERLAPMATRIX Area-weighting matrix for 1D resampling
% Returns a sparse nOut x nIn matrix M, where M(j,k) is the length of the
% overlap between output pixel j and input pixel k divided by psOut. Both
% grids are centered at the same point. The rows sum to the fraction of the
% output pixel covered by the input grid (1 inside, less than 1 at the
% edges, 0 outside).
nIn = double(nIn);
nOut = double(nOut);
edgesIn = (2 * (0:nIn) - nIn) * psIn / 2;
edgesOut = (2 * (0:nOut) - nOut) * psOut / 2;
ii = cell(nOut,1);
jj = cell(nOut,1);
vv = cell(nOut,1);
for j = 1:nOut
    % Candidate input pixels (one extra on both sides, the overlap is
    % clipped to zero below)
    kLo = max(1, floor((edgesOut(j) - edgesIn(1)) / psIn));
    kHi = min(nIn, ceil((edgesOut(j+1) - edgesIn(1)) / psIn) + 1);
    k = kLo:kHi;
    ov = max(0, min(edgesOut(j+1), edgesIn(k+1)) - max(edgesOut(j), edgesIn(k)));
    ii{j} = j * ones(1, numel(k));
    jj{j} = k;
    vv{j} = ov / psOut;
end
M = sparse([ii{:}], [jj{:}], [vv{:}], nOut, nIn);