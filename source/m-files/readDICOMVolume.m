function [V, spatial] = readDICOMVolume(fpath)
%READDICOMVOLUME Read a DICOM series from a folder as a volume
%   [V, SPATIAL] = READDICOMVOLUME(FPATH) reads all DICOM files of the
%   folder FPATH. If dicomreadVolume is available it is used. Otherwise a
%   fallback based on dicominfo and dicomread is used (e.g. GNU Octave).
%   The slices are sorted along the slice normal (from
%   ImageOrientationPatient, or by the z position if that is missing).
%
%   V is of size [Rows Columns 1 nSlices]. SPATIAL has the fields
%   PatientPositions (nSlices x 3), PixelSpacings (nSlices x 2, [row column]
%   spacing), PatientOrientations (2 x 3 x nSlices, only if
%   ImageOrientationPatient exists) and ImageSize ([Rows Columns nSlices]).
%
% See also dicomreadVolume, dicominfo, dicomread

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2026 Ville-Veikko Wettenhovi
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('dicomreadVolume', 'file') == 2
    [V, spatial] = dicomreadVolume(fpath);
    return
end

files = dicomFiles(fpath);
nSlices = numel(files);
if nSlices == 0
    error('No DICOM files found in %s!', fpath);
end

pos = zeros(nSlices, 3);
spacing = zeros(nSlices, 2);
orient = zeros(2, 3, nSlices);
hasIOP = true;
imgs = cell(nSlices, 1);
for kk = 1 : nSlices
    f = fullfile(fpath, files(kk).name);
    try
        info = dicominfo(f);
    catch ME
        error('Could not read the DICOM file %s: %s', f, ME.message);
    end
    imgs{kk} = dicomread(info);
    if ~isfield(info, 'ImagePositionPatient') || numel(info.ImagePositionPatient) < 3
        error('The DICOM file %s does not contain ImagePositionPatient!', f);
    end
    if ~isfield(info, 'PixelSpacing') || numel(info.PixelSpacing) < 2
        error('The DICOM file %s does not contain PixelSpacing!', f);
    end
    pos(kk, :) = double(info.ImagePositionPatient(1 : 3))';
    spacing(kk, :) = double(info.PixelSpacing(1 : 2))';
    if hasIOP && isfield(info, 'ImageOrientationPatient') && numel(info.ImageOrientationPatient) >= 6
        iop = double(info.ImageOrientationPatient(1 : 6));
        orient(:, :, kk) = reshape(iop, 3, 2)';
    else
        hasIOP = false;
    end
    if kk > 1 && ~isequal(size(imgs{kk}), size(imgs{1}))
        error('The slice sizes of the DICOM files in %s differ!', fpath);
    end
end

% Sort along the slice normal
if hasIOP
    normal = cross(orient(1, :, 1), orient(2, :, 1));
    dist = pos * normal(:);
else
    dist = pos(:, 3);
end
[~, order] = sort(dist);

nRows = size(imgs{1}, 1);
nCols = size(imgs{1}, 2);
V = zeros(nRows, nCols, 1, nSlices, class(imgs{1}));
for kk = 1 : nSlices
    V(:, :, 1, kk) = imgs{order(kk)};
end

spatial.PatientPositions = pos(order, :);
spatial.PixelSpacings = spacing(order, :);
if hasIOP
    spatial.PatientOrientations = orient(:, :, order);
end
spatial.ImageSize = [nRows nCols nSlices];
end

function files = dicomFiles(fpath)
% DICOM files of the folder: all files with the extension .dcm or .ima
% (case-insensitive) and all other files (e.g. without an extension) that
% isdicom accepts, i.e. the files that dicomreadVolume reads
    files = dir(fpath);
    files = files(~[files.isdir]);
    [~, ~, ext] = cellfun(@fileparts, {files.name}, 'UniformOutput', false);
    known = ismember(lower(ext), {'.dcm', '.ima'});
    for ii = find(~known(:))'
        known(ii) = isdicom(fullfile(fpath, files(ii).name));
    end
    files = files(known);
end
