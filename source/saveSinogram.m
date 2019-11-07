function saveSinogram(sino, varargin)
%SAVESINOGRAM Saves the specified sinogram in the specified format
%   This function saves the input sinogram to the specified format.
%   Available formats are DICOM, NIfTI, Analyze 7.5 or raw binary image.
%   DICOM support requires image processing toolbox (or dicom package on
%   Octave (untested)), NIfTI and Analyze support require "Tools for NIfTI
%   and ANALYZE image" toolbox from MathWorks file exchange: 
%   https://se.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
%
%   Raw data format is used by default. DICOM format does not support
%   time-series of images.
%
% Examples:
%   saveSinogram(sino)
%   saveSinogram(sino, 'nifti', 'outputFileName', options)
%
% Inputs:
%   sino = 3D, 4D or a cell matrix containing the sinograms. This can be
%   e.g. the SinM-matrix output by form_sinograms.m
%
%   If a 4D matrix is used, the fourth dimension is assumed to be the
%   number time steps.
%
%   If a cell matrix is used, then each separate cell is assumed to be
%   another time step.
%
%   type = The output file type. 'nifti' uses NIfTI, 'analyze' uses Analyze
%   7.5 format, 'dicom' uses DICOM format and 'raw' uses raw binary data
%   (default).
%
%   Filename = The filename (and optionally full path) of the output file.
%   If omitted or an empty array is used, the same naming scheme will be
%   used as elsewhere in OMEGA. Omitting this will require the usage of the
%   below options struct.
%
%   options = The struct created by any of the main-files, e.g.
%   gate_main.m. Used for naming purposes and to save the sinogram bin
%   sizes.
%
%   Using DICOM format creates a single file for each slice.
%
% Outputs:
%   The output file will be saved either in the current directory (if no
%   filename was specified) or in the path specified in the filename.
%
% See also

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

if nargin >= 3 && ~isempty(varargin{2})
    filename = varargin{1};
elseif nargin >= 4 && ~isempty(varargin{3})
    filename = [varargin{3}.machine_name '_' varargin{3}.name];
else
    error('No output filename specified')
end
if iscell(sino)
    n_time_steps = length(sino);
    sino = reshape(sino, 1, 1, 1, n_time_steps);
    sino = cell2mat(sino);
end
if nargin >= 2 && ~isempty(varargin{1})
    type = varargin{1};
else
    type = 'raw';
end

origin = [0 0 0];
description = 'OMEGA sinogram';
kuva = sino;
if size(kuva,2) == 1
    error('The size of the input sinogram has to be at least XxYxZ')
end
if nargin >= 4 && ~isempty(varargin{3})
    voxel_size = [varargin{3}.cr_p/2 varargin{3}.cr_p/2 varargin{3}.TotSinos/size(kuva,3)];
else
    voxel_size = [0 0 0];
end
if strcmp(type, 'nifti')
    filename = [filename '.nii'];
    kuva = rot90(kuva,1);
    if exist('niftiwrite', 'file') == 2
        niftiwrite(kuva,filename);
        info = niftiinfo(filename);
        info.Description = description;
        info.PixelDimensions = voxel_size;
        niftiwrite(kuva,filename,info);
    else
        nii = make_nii(kuva, voxel_size, origin, [], description);
        save_nii(nii, filename);
    end
elseif strcmp(type, 'analyze')
    filename = [filename '.img'];
    ana = make_ana(kuva, voxel_size, origin, [], description);
    save_untouch_nii(ana, filename);
elseif strcmp(type, 'dicom')
    if size(kuva,4) > 1
        error('4D data is not supported by DICOM conversion')
    end
    kuva = double(kuva);
    for ll = 1 : prop.Nz
        dicomwrite(kuva(:,:,ll), [filename '_slice' num2str(ll) '.dcm']);
    end
elseif strcmp(type, 'interfile')
    saveInterfile(filename, kuva, 'sinogram', class(kuva));
elseif strcmp(type, 'metaimage')
    saveMetaImage(filename, kuva);
elseif strcmp(type, 'raw')
    filename = [filename '.raw'];
    fid = fopen(filename);
    fwrite(fid, kuva(:),class(kuva));
    fclose(fid);
else
    error('Unsupported output filetype')
end

