function options = loadInveonCTData(options,varargin)
%LOADINVEONCTDATA Loads Inveon CT projection data
%   Loads the Siemens Inveon CT CAT-data files by using the corresponding
%   header file. The user will be prompted for the input header file.
%   Flat-field (dark/light normalization) correction will be performed
%   automatically. Binning can also be included. The final images have been
%   flipped in the horizontal direction.
%   
%   Examples:
%       options = loadInveonCTData(options);
%       options = loadInveonCTData(options, saveData, fileName);
%   Inputs:
%       options = All the required variables will be saved to the
%       options-struct, including the projection images.
%       saveData = (Optional) If set to true, will save the projection
%       images in a mat file with the filename specified with fileName.
%       fileName = (Optional) The name of the mat-file where the projection
%       images are stored. Required if saveData is set to true.
%
% See also loadProjectionData, loadProjectionImages
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021 Ville-Veikko Wettenhovi
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

if nargin >= 2 && ~isempty(varargin) && ~isempty(varargin{1}) && ~isempty(varargin{2})
    saveData = varargin{1};
    fileName = varargin{2};
else
    saveData = false;
end
[file, fpath] = uigetfile('*.hdr','Select Inveon CT projection (CAT) header file');
if isequal(file, 0)
    error('No file was selected')
end
nimi = [fpath file];
fid = fopen(nimi);
hdr = textscan(fid,'%s','CommentStyle','#');
hdr = hdr{1};
fclose(fid);
% hdr = textread(nimi, '%s','commentstyle','shell');
options.nBed = str2double(hdr{find(strcmp(hdr,'number_of_bed_positions')) + 1});
if numel(options.use_N_positions) > options.nBed || max(options.use_N_positions) > options.nBed
    error(['The number of selected bed positions is higher than the number of bed positions measured! '...
        num2str(options.nBed) ' bed positions available.'])
end
% trBinSize = str2double(hdr{find(strcmp(hdr,'transaxial_bin_size')) + 1});
% axPlaneSize = str2double(hdr{find(strcmp(hdr,'axial_plane_size')) + 1});
dType = str2double(hdr{find(strcmp(hdr,'data_type')) + 1});
options.nRowsD = str2double(hdr{find(strcmp(hdr,'x_dimension')) + 1});
options.nColsD = str2double(hdr{find(strcmp(hdr,'y_dimension')) + 1});
options.zSize = str2double(hdr{find(strcmp(hdr,'z_dimension')) + 1});
angleStart = str2double(hdr{find(strcmp(hdr,'rotating_stage_start_position')) + 1});
angleStop = str2double(hdr{find(strcmp(hdr,'rotating_stage_stop_position')) + 1});
options.nProjections = str2double(hdr{find(strcmp(hdr,'number_of_projections')) + 1});
% uncroppedTrPixels = str2double(hdr{find(strcmp(hdr,'ct_uncropped_transaxial_pixels')) + 1});
% uncroppedAxPixels = str2double(hdr{find(strcmp(hdr,'ct_uncropped_axial_pixels')) + 1});
% croppedTrPixels = str2double(hdr{find(strcmp(hdr,'ct_cropped_transaxial_pixels')) + 1});
% croppedAxPixels = str2double(hdr{find(strcmp(hdr,'ct_cropped_axial_pixels')) + 1});
options.dPitch = str2double(hdr{find(strcmp(hdr,'ct_xray_detector_pitch')) + 1})/10^3;
options.sourceToDetector = str2double(hdr{find(strcmp(hdr,'ct_source_to_detector')) + 1});
options.sourceToCRot = str2double(hdr{find(strcmp(hdr,'ct_source_to_crot')) + 1});
bedOffset = str2double(hdr(find(strcmp(hdr,'ending_bed_offset')) + 1));
options.bedOffset = bedOffset(:) * 10;
if numel(options.use_N_positions) < options.nBed
    options.bedOffset = options.bedOffset(options.use_N_positions);
end
averageOffset = str2double(hdr(find(strcmp(hdr,'ct_projection_average_center_offset')) + 1));
options.bedOffset(2:end) = options.bedOffset(2:end) - averageOffset * (1 : numel(options.use_N_positions) - 1)';
options.bedOffset = options.bedOffset - options.bedOffset(1);
% verticalOffset = str2double(hdr(find(strcmp(hdr,'ct_projection_average_center_offset')) + 1));
horizontalOffset = str2double(hdr(find(strcmp(hdr,'ct_projection_center_offset')) + 2));
horizontalBedOffset = str2double(hdr(find(strcmp(hdr,'ct_projection_horizontal_bed_offset')) + 2));
options.horizontalOffset = horizontalOffset - horizontalBedOffset;
% binningTr = croppedTrPixels / nColsD;
% binningAx = croppedAxPixels / nColsD;

options.angles = -linspace(angleStart,angleStop,options.nProjections)';

fid = fopen([fpath file(1:end-4)]);
if dType == 2
    tyyppi = 'int16=>int16';
elseif dType == 3
    tyyppi = 'int32=>int32';
elseif dType == 4
    tyyppi = 'single=>single';
else
    error('Unknown file type')
end
projData = fread(fid, inf, tyyppi);
fclose(fid);
C = projData((numel(projData) - options.nColsD*options.nRowsD*options.zSize*options.nBed - options.nColsD*options.nRowsD*2) + 1:...
    numel(projData) - options.nColsD*options.nRowsD*options.zSize*options.nBed);
dark = C(1:end/2);
light = C(end/2+1:end) - dark;
dark = reshape(dark, options.nRowsD, options.nColsD);
light = reshape(light, options.nRowsD, options.nColsD);
projData = projData((numel(projData) - options.nColsD*options.nRowsD*options.zSize*options.nBed) + 1:end);
projData = reshape(projData,options.nRowsD,options.nColsD,options.zSize,options.nBed);
if numel(options.use_N_positions) ~= options.nBed
    projData = projData(:,:,:,options.use_N_positions);
end
options.nProjections = options.nProjections * numel(options.use_N_positions);
projData = bsxfun(@minus, projData, dark);
projData = bsxfun(@rdivide,single(projData),single(light));
if options.binning > 1
    B = sum(reshape(projData,options.binning,[]),1,'native');
    B = squeeze(sum(reshape(B,size(projData,1)/options.binning,options.binning,[]),2,'native'));
    B = reshape(B, size(projData,1)/options.binning,size(projData,2)/options.binning,options.nProjections,numel(options.use_N_positions));
    projData = B;
    options.nRowsD = size(projData,1);
    options.nColsD = size(projData,2);
end
projData = reshape(projData, size(projData,1), size(projData,2), size(projData,3) * numel(options.use_N_positions));
projData = fliplr(projData);
if saveData
    if exist('OCTAVE_VERSION','builtin') == 0
        save(fileName,'projData','-v7.3')
    else
        save(fileName,'projData','-v7')
    end
end
if numel(options.use_N_positions) < options.nBed
    options.nBed = numel(options.use_N_positions);
end
options.SinM = projData;

end

