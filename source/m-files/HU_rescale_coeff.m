function [HU_a, HU_b] = HU_rescale_coeff(fpathCT, patientPositions)
%HU_RESCALE_COEFF DICOM rescale slope and intercept of each slice of a CT series
%   [HU_a, HU_b] = HU_RESCALE_COEFF(fpathCT) returns the RescaleSlope and
%   RescaleIntercept of the .dcm files in the folder fpathCT in file name
%   order, as 1 x 1 x nFiles arrays.
%
%   [HU_a, HU_b] = HU_RESCALE_COEFF(fpathCT, patientPositions) returns the
%   coefficients in the slice order of dicomreadVolume. patientPositions is
%   spatial.PatientPositions (nSlices x 3) from dicomreadVolume. Each file
%   is matched to the slice whose position is closest to the file's
%   ImagePositionPatient. This is required when the file names are not in
%   slice order. Files with the extension .dcm or .ima (case-insensitive)
%   and other DICOM files (e.g. without an extension) are used. An error is thrown if the number of files differs from the
%   number of slices or if the files do not map one-to-one to the slices.
    if nargin < 2 || isempty(patientPositions)
        files = dicomFiles(fpathCT);
        n = numel(files);
        HU_a = zeros(n, 1);
        HU_b = zeros(n, 1);
        for ii = 1:n
            f = fullfile(fpathCT, files(ii).name);
            info = dicominfo(f);
            HU_a(ii) = double(info.RescaleSlope);
            HU_b(ii) = double(info.RescaleIntercept);
        end
        HU_a = reshape(HU_a, 1, 1, n);
        HU_b = reshape(HU_b, 1, 1, n);
        return
    end
    files = dicomFiles(fpathCT);
    n = numel(files);
    nSlices = size(patientPositions, 1);
    if n ~= nSlices
        error('The number of DICOM files (.dcm or .ima) in %s (%d) differs from the number of CT slices (%d)!', fpathCT, n, nSlices);
    end
    slope = zeros(n, 1);
    intercept = zeros(n, 1);
    sliceIdx = zeros(n, 1);
    for ii = 1:n
        info = dicominfo(fullfile(fpathCT, files(ii).name));
        slope(ii) = double(info.RescaleSlope);
        intercept(ii) = double(info.RescaleIntercept);
        pos = double(info.ImagePositionPatient(:))';
        [~, sliceIdx(ii)] = min(sum((patientPositions - pos).^2, 2));
    end
    if numel(unique(sliceIdx)) ~= n
        error('The DICOM files in %s do not map one-to-one to the CT slices (two files have the same position)!', fpathCT);
    end
    HU_a = zeros(n, 1);
    HU_b = zeros(n, 1);
    HU_a(sliceIdx) = slope;
    HU_b(sliceIdx) = intercept;
    HU_a = reshape(HU_a, 1, 1, n);
    HU_b = reshape(HU_b, 1, 1, n);
end

function files = dicomFiles(fpathCT)
% DICOM files of the folder: all files with the extension .dcm or .ima
% (case-insensitive) and all other files (e.g. without an extension) that
% isdicom accepts, i.e. the files that dicomreadVolume reads
    files = dir(fpathCT);
    files = files(~[files.isdir]);
    [~, ~, ext] = cellfun(@fileparts, {files.name}, 'UniformOutput', false);
    known = ismember(lower(ext), {'.dcm', '.ima'});
    for ii = find(~known(:))'
        known(ii) = isdicom(fullfile(fpathCT, files(ii).name));
    end
    files = files(known);
end
