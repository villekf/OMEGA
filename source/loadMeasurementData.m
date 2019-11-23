function options = loadMeasurementData(options, varargin)
%LOADMEASUREMENTDATA Loads the measurement data
%   Loads the measurement data for non-GATE situations from either a
%   mat-file, a NIfTI file, Analyze 7.5 file, DICOM file, Interfile or raw
%   data file 
%   in the specified bit-format. For the last case the default is 32-bit
%   integer (int32). A different bit format for the raw data file can be
%   input as an optional parameter. The measurement data will be saved in
%   structure as SinM for sinogram data or coincidences for raw list-mode
%   data. Raw list-mode data supports only mat-files.
%
%   Loading DICOM files requires image processing toolbox on MATLAB or the
%   dicom-package on Octave (untested). For NIfTI and Analyze format the
%   "Tools for NIfTI and ANALYZE image" toolbox from MathWorks file
%   exchange is needed
%   (https://se.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
%
%   Measurement data previously loaded (i.e. either SinM or coincidences
%   exist in the options struct) can be used if no data is selected when
%   prompted.
%
%   The file extension is used to determine the input data format. E.g.
%   img/hdr/nii are considered Analyze or NIfTI format files. For
%   Interfile i33/h33 are considered Interfile files.
%
% Examples:
%   options = loadMeasurementData(options)
%   options = loadMeasurementData(options,'uint16', skip)
%
% Inputs:
%   type = The type of raw data used. Used only for raw data. E.g. 'int8',
%   'uint16', 'int32', 'single', 'double', 'uint64'.
%
%   skip = Number of bytes skipped before the actual data. I.e. the length
%   of possible header in bytes. If the header is at the end of file, use
%   negative value. Unnecessary if no header is present. Applicable ONLY
%   for raw data.
%
% See also dicomread, load_nii, niftiread, loadInterfile
%

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
[file, fpath] = uigetfile('*.*','Select PET Sinogram or list-mode data');

if isequal(file, 0)
    if isfield(options,'SinM') || isfield(options,'coincidences')
        warning('Previously loaded measurement data found. Using the data already stored in struct.')
        return
    else
        error('No file was selected')
    end
end

FileName = fullfile(fpath, file);

if nargin > 1
    type = varargin{1};
else
    type = 'int32';
end

if nargin > 2
    skip = varargin{2};
else
    skip = 0;
end

if strcmp(FileName(end-2:end),'mat')
    storedStructure = load(FileName);
    variables = fields(storedStructure);
    
    if options.use_raw_data == false
        if isfield(storedStructure, 'SinM')
            options.SinM = storedStructure.SinM;
        elseif isfield(storedStructure, 'raw_SinM')
            options.SinM = storedStructure.raw_SinM;
        else
            options.SinM = storedStructure.(variables{1});
            if length(variables) > 1
                warning('The input data file contains more than one variable! The first one is used')
            end
        end
    else
        if isfield(storedStructure, 'coincidences')
            options.coincidences = storedStructure.coincidences;
        else
            options.coincidences = storedStructure.(variables{1});
            if length(variables) > 1
                warning('The input data file contains more than one variable! The first one is used')
            end
        end
    end
elseif strcmp(FileName(end-2:end),'img') || strcmp(FileName(end-2:end),'nii') || strcmp(FileName(end-2:end),'hdr') || strcmp(FileName(end-2:end),'.gz')
    if exist('niftiread','file') == 2
        img = niftiread(FileName);
    elseif exist('load_nii', 'file') == 2
        nii = load_nii(FileName);
        img = nii.img;
    else
        error('Image processing toolbox or Tools for NIfTI and ANALYZE image toolbox not found, check that the files are in MATLAB path')
    end
    if size(img,1) ~= options.Ndist
        if size(img,2) ~= options.Nang
            error('Size mismatch between input file and sinogram dimensions')
        else
            img = flipud(permute(img, [2 1 3 4]));
        end
    end
    if options.partions > 1
        if size(img,4) > 1
            options.SinM = cell(options.partitions,1);
            for kk = 1 : options.partitions
                options.SinM{kk} = img(:,:,:,kk);
            end
        else
            error('No 4D data found')
        end
    else
        options.SinM = zeros(options.Ndist, options.Nang, options.TotSinos, 'single');
        for kk = 1 : size(img,4)
            options.SinM = options.SinM + img(:,:,:,kk);
        end
    end
elseif strcmp(FileName(end-2:end),'dcm') || strcmp(FileName(end-2:end),'dicom')
    if exist('dicomread', 'file') == 2
        img = squeeze(dicomread(FileName));
        if size(img,1) ~= options.Ndist
            if size(img,2) ~= options.Nang
                error('Size mismatch between input file and sinogram dimensions')
            else
                img = permute(img, [2 1 3 4]);
            end
        end
        if options.partions > 1
            error('No 4D data in DICOM files')
        end
        options.SinM = img;
    else
        error('No dicomread file found. Image processing toolbox is required for reading DICOM files on MATLAB')
    end
elseif strcmp(FileName(end-2:end),'i33') || strcmp(FileName(end-2:end),'h33')
    options.SinM = loadInterfile(FileName);
elseif strcmp(FileName(end-2:end),'mhd') || strcmp(FileName(end-2:end),'mha')
    options.SinM = loadMetaImage(FileName);
else
    fid = fopen(nimi);
    f_size = dir(nimi);
    f_size = f_size.bytes;
    if nargin >= 2 && skip > 0
        fread(fid,skip,'*uint8');
        f_size = f_size - skip;
    elseif nargin >= 2 && skip < 0
        f_size = f_size + skip;
    end
    if strcmp(type,'int16') || strcmp(type,'uint16')
        f_size = f_size / 2;
    elseif strcmp(type,'int32') || strcmp(type,'uint32') || strcmp(type,'single')
        f_size = f_size / 4;
    elseif strcmp(type,'int64') || strcmp(type,'uint64') || strcmp(type,'double')
        f_size = f_size / 8;
    end
    koko = options.Nang*options.Ndist*options.TotSinos;
    if (f_size > koko || f_size < koko) && options.partitions == 1
        if mod(f_size, koko) == 0 && f_size > koko
            warning('4th dimension will be merged')
        else
            fclose(fid);
            error(['Size mismatch between input file and sinogram dimensions. Expected [' num2str(options.Ndist) ' ' ...
                num2str(options.Nang) ' ' num2str(options.TotSinos) '].'])
        end
    elseif (f_size > koko * options.partitions || f_size < koko * options.partitions)
        fclose(fid);
        error(['Size mismatch between input file and sinogram dimensions. Expected [' num2str(options.Ndist) ' ' ...
            num2str(options.Nang) ' ' num2str(options.TotSinos) ' ' num2str(options.partitions) '].'])
    end
    img = fread(fid, [options.Ndist, options.Nang, options.TotSinos, f_size/koko], [type '=>single'],0,'l');
    if options.partions > 1
        if size(img,4) > 1
            options.SinM = cell(options.partitions,1);
            for kk = 1 : options.partitions
                options.SinM{kk} = img(:,:,:,kk);
            end
        else
            error('No 4D data found')
        end
    else
        options.SinM = zeros(options.Ndist, options.Nang, options.TotSinos, 'single');
        for kk = 1 : size(img,4)
            options.SinM = options.SinM + img(:,:,:,kk);
        end
    end
end
end

