function saveImage(img, varargin)
%SAVEIMAGE Saves the specified reconstructed images in the specified format
%   This function saves the input image(s) to the specified format.
%   Available formats are DICOM, NIfTI, Analyze 7.5, Interfile, MetaImage
%   or raw image.
%
%   DICOM support requires image processing toolbox (or dicom
%   package on Octave (untested)), NIfTI requires either image processing
%   toolbox (R2017b and up) OR "Tools for NIfTI and ANALYZE image" toolbox
%   from MathWorks file exchange:
%   https://se.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
%   Analyze format requires "Tools for NIfTI and ANALYZE image". MetaImage
%   requires "Medical Image Processing Toolbox" from file exchange:
%   https://se.mathworks.com/matlabcentral/fileexchange/41594-medical-image-processing-toolbox
%
%   Raw data format is used by default. DICOM format does not support
%   time-series of images. Interfile, MetaImage and raw format are the
%   same, except Interfile and MetaImage also output the header file.
%
% Examples:
%   saveImage(img)
%   saveImage(img, type, filename, options)
%   saveImage(img, 'nifti', 'outputFileName', options)
%
% Inputs:
%   img = 3D, 4D, 5D or a cell matrix containing the estimates. This can be
%   e.g. the pz-matrix output by reconstructions_main
%
%   If a 4D matrix is used, the fourth dimension is assumed to be the
%   number of iterations. The last iteration is used in this case.
%
%   If a 5D matrix is used, the first four dimensions behave as with 4D
%   matrix, but the 5th dimension is considered as time.
%
%   If a cell matrix is used, then all non-empty cells are saved in
%   separate files. The last iterations will be used and in the case of
%   time-series data, the time dimension will also be saved (if the output
%   format supports it).
%
%   type = The output file type. 'nifti' uses NIfTI, 'analyze' uses Analyze
%   7.5 format, 'dicom' uses DICOM format, 'interfile' uses Interfile
%   format, 'metaimage' uses MetaImage, and 'raw' uses raw binary data.
%
%   DICOM uses double precision regardless of the input image. Other
%   formats use the same type as the input image.
%
%   Using DICOM format creates a single file for each slice.
%
%   Filename = The filename (and optionally full path) of the output file.
%   If omitted or an empty array is used, the same naming scheme will be
%   used as elsewhere in OMEGA. Omitting this will require the usage of a
%   cell-matrix or the below options struct to be present.
%
%   options = The struct created by any of the main-files, e.g. gate_main.m
%   Omitting this and not using the cell matrix created by OMEGA prevents
%   the saving of additional information in the output files (e.g. voxel
%   size)
%
% Outputs:
%   The output file will be saved either in the current directory (if no
%   filename was specified) or in the path specified in the filename.
%
% See also saveInterfile, importData, saveSinogram, niftiwrite, save_nii

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ville-Veikko Wettenhovi
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

if iscell(img)
    prop = img{end,1};
end
if nargin >= 3 && ~isempty(varargin{2})
    orig_filename = varargin{2};
else
    if iscell(img)
        orig_filename = [prop.machine_name '_' prop.name];
    else
        if nargin >= 4 && ~isempty(varargin{3})
            orig_filename = [varargin{3}.machine_name '_' varargin{3}.name];
        else
            error('No output filename found. Use "saveImage(img, ''outputFileName'')" to set the output filename')
        end
    end
end
if nargin >= 2 && ~isempty(varargin{1})
    type = varargin{1};
else
    type = 'raw';
end

if iscell(img)
    pituus = size(img,1) - 1;
else
    pituus = 1;
end
% Loop through all the algorithms
for kk = 1 : pituus
    origin = [0 0 0];
    description = 'OMEGA reconstruction';
    if iscell(img)
        if isempty(img{kk,1})
            continue
        end
        % Dynamic case
        if prop.n_time_steps > 1
            kuva = img(kk,:);
            kuva = reshape(kuva, 1, 1, 1, 1, prop.n_time_steps);
            kuva = cell2mat(kuva);
        else
            kuva = img{kk};
        end
        % Last iteration
        kuva = squeeze(kuva(:,:,:,end,:));
        voxel_size = [prop.FOV_x/prop.Nx prop.FOV_y/prop.Ny prop.axial_FOV/prop.Ny];
        reko = algorithms_char();
        reko = reko{kk};
        description = [description ' with ' reko];
        filename = [orig_filename '_' reko];
    else
        kuva = img;
        % Last iteration
        if size(kuva,5) > 1
            kuva = squeeze(kuva(:,:,:,end,:));
        elseif size(kuva,4) > 1
            kuva = squeeze(kuva(:,:,:,end,:));
        elseif size(kuva,2) == 1
            error('The size of the input image has to be at least XxYxZ')
        end
        if nargin >= 4 && ~isempty(varargin{3})
            voxel_size = [varargin{3}.FOVa_x/varargin{3}.Nx varargin{3}.FOVa_y/varargin{3}.Ny varargin{3}.axial_fov/varargin{3}.Ny];
        else
            voxel_size = [0 0 0];
        end
        if nargin >= 4 && ~isempty(varargin{3})
            if isfield(varargin{3}, 'simple')
                description = [description ' with OSEM'];
                filename = [orig_filename '_OSEM'];
            else
                filename = orig_filename;
            end
        else
            filename = orig_filename;
            reko = [];
        end
    end
    % NIfTI
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
        % Analyze 7.5
    elseif strcmp(type, 'analyze')
        filename = [filename '.img'];
        ana = make_ana(kuva, voxel_size, origin, [], description);
        save_untouch_nii(ana, filename);
        % DICOM
    elseif strcmp(type, 'dicom')
        if size(kuva,4) > 1
            error('4D data is not supported by DICOM conversion')
        end
        kuva = double(kuva);
        for ll = 1 : size(kuva,3)
            dicomwrite(kuva(:,:,ll), [filename '_slice' num2str(ll) '.dcm']);
        end
        % Interfile
    elseif strcmp(type,'interfile')
        if iscell(img)
            saveInterfile(filename, kuva, reko, class(kuva), prop)
        elseif nargin >= 4 && ~isempty(varargin{3})
            saveInterfile(filename, kuva, reko, class(kuva), varargin{3})
        else
            saveInterfile(filename, kuva, reko, class(kuva))
        end
        % MetaImage
    elseif strcmp(type, 'metaimage')
        saveMetaImage(filename,kuva);
        % Binary
    elseif strcmp(type, 'raw')
        filename = [filename '.raw'];
        fid = fopen(filename);
        fwrite(fid, kuva(:),class(kuva));
        fclose(fid);
    else
        error('Unsupported output filetype')
    end
end

