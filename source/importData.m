function output = importData(type, varargin)
%IMPORTDATA Import data into MATLAB matrix
%   Import the input data (user will be prompted for the data) of the
%   specified type into a MATLAB matrix.
%
%   Available types are DICOM, NIfTI, Analyze 7.5, Interfile, MetaImage and
%   raw data format.
%
%   Loading DICOM files requires image processing toolbox on MATLAB or the
%   dicom-package on Octave (untested). NIfTI requires either image
%   processing toolbox or the "Tools for NIfTI and ANALYZE image" toolbox
%   from MathWorks file exchange. Analyze requires "Tools for NIfTI and
%   ANALYZE image" toolbox
%   (https://se.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
%
% Examples:
%   output = importData('nifti')
%   output = importData('raw','uint32',N1,N2,N3,N4,N5,skip)
%
% Inputs:
%   type = The type of data used. Available ones are 'dicom', 'nifti',
%   'analyze', 'interfile', 'metaimage' and 'raw'.
%
%   format = The data format type of the raw data. E.g. 'single', 'uint32',
%   'double', etc.
%
%   N1,...,N5 = The size of the specific dimensions. E.g. if N3 = 64, then
%   64 slices are assumed to be in the third dimension.. Up to five
%   dimensions are supported. This is ONLY used with raw data. On other
%   formats, the dimensions are extracted from the headers/metainfo. If a
%   specific dimension does not exist, use 1 or [].
%
%   skip = Number of bytes skipped before the actual data. I.e. the length
%   of possible header in bytes. If the header is at the end of file, use
%   negative value. Unnecessary if no header is present. Applicable ONLY
%   for raw data. This is computed automatically if the header is in the
%   beginning of the file and the data format and dimensions are input.
%
% See also dicomread, niftiread, load_nii, loadInterfile, analyze75read
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

if strcmp(type,'nifti')
    [file, fpath] = uigetfile(('*.nii;*.img;*.hdr;*.gz'),'Select NIfTI data file');
    FileName = fullfile(fpath, file);
    if exist('niftiread','file') == 2
        output = niftiread(FileName);
    elseif exist('load_nii', 'file') == 2
        nii = load_nii(FileName);
        output = nii.img;
    else
        error('Image processing toolbox or Tools for NIfTI and ANALYZE image toolbox not found, check that the files are in MATLAB path')
    end
    
elseif strcmp(type, 'dicom')
    [file, fpath] = uigetfile(('*.dcm'),'Select DICOM data file');
    FileName = fullfile(fpath, file);
    if exist('dicomread', 'file') == 2
        output = squeeze(dicomread(FileName));
    else
        error('No dicomread file found. Image processing toolbox is required for reading DICOM files on MATLAB')
    end
elseif strcmp(type, 'analyze')
    [file, fpath] = uigetfile(('*.img;*.hdr'),'Select Analyze data file');
    FileName = fullfile(fpath, file);
    if exist('analyze75read','file') == 2
        output = analyze75read(FileName);
    elseif exist('load_nii', 'file') == 2
        nii = load_nii(FileName);
        output = nii.img;
    else
        error('Image processing toolbox or Tools for NIfTI and ANALYZE image toolbox not found, check that the files are in MATLAB path')
    end
    
elseif strcmp(type, 'interfile')
    [file, fpath] = uigetfile(('*.img;*.hdr;*.i33;*.h33'),'Select Interfile data file');
    FileName = fullfile(fpath, file);
    output = loadInterfile(FileName);
elseif strcmp(type, 'metaimage')
    [file, fpath] = uigetfile(('*.mhd;*.mha;*.raw'),'Select MetaImage header file');
    FileName = fullfile(fpath, file);
    output = loadMetaImage(FileName);
elseif strcmp(type, 'raw')
    type = varargin{1};
    n_bytes = struct('double',8,'single', 4, 'int8',1,'uint8',1,'int16',2,'uint16',2,'int32',4,'uint32',4, 'int64',8, 'uint64', 8);
    if nargin >= 3 && ~isempty(varargin{2})
        N1 = varargin{2};
    else
        N1 = 1;
    end
    if nargin >= 4 && ~isempty(varargin{3})
        N2 = varargin{3};
    else
        N2 = 1;
    end
    if nargin >= 5 && ~isempty(varargin{4})
        N3 = varargin{4};
    else
        N3 = 1;
    end
    if nargin >= 6 && ~isempty(varargin{5})
        N4 = varargin{5};
    else
        N4 = 1;
    end
    if nargin >= 7 && ~isempty(varargin{6})
        N5 = varargin{6};
    else
        N5 = 1;
    end
    if nargin >= 8 && ~isempty(varargin{7})
        skip = varargin{7};
    else
        skip = 0;
    end
    [file, fpath] = uigetfile('*.*','Select raw data file');
    FileName = fullfile(fpath, file);
    f_size = dir(FileName);
    f_size = f_size.bytes;
    fid = fopen(FileName);
    if nargin >= 8 && skip > 0
        fread(fid,skip,'*uint8');
        f_size = f_size - skip;
    elseif nargin >= 8 && skip < 0
        f_size = f_size + skip;
    else
        skip = f_size - n_bytes.(type) * N1 * N2 * N3 * N4 * N5;
        if skip > 0
            fread(fid,skip,'*uint8');
            f_size = f_size - skip;
        end
    end
    if strcmp(type,'int16') || strcmp(type,'uint16')
        f_size = f_size / 2;
    elseif strcmp(type,'int32') || strcmp(type,'uint32') || strcmp(type,'single')
        f_size = f_size / 4;
    elseif strcmp(type,'int64') || strcmp(type,'uint64') || strcmp(type,'double')
        f_size = f_size / 8;
    end
    output = fread(fid, f_size, [type '=>' type]);
    fclose(fid);
    output = reshape(output, N1, N2, N3, N4, N5);
end

