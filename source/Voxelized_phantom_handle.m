function A = Voxelized_phantom_handle(input_name, output_name, pixdim, pixsize, varargin)
% VOXELIZED_PHANTOM_HANDLE This code can be used to crop a voxelized
% phantom and automatically convert it to uint16 format in Interfile or
% MetaImage file format.
%
%   A voxelized phantom can be e.g. XCAT phantom or a collection of DICOM
%   CT-images. This function automatically crops the phantoms such that
%   there are no empty rows/columns/slices. Supports both static and
%   dynamic data, though with dynamic data there are a few air gaps in the
%   phantom to make sure there is no cut-off due to motion. All values are 
%   multiplied such that the maximum value will be 1000 and then converted
%   to unsigned 16-bit integer format. These values can then be used to
%   form the Attrange.dat file. The phantom is saved in Interfile or
%   MetaImage format with automatically created header file that can be
%   input into GATE. Negative values are removed by adding the smallest
%   (largest negative) value to the phantom. The header files use Unix file
%   format by default, but this can be (optionally) altered to Windows file
%   format instead. Interfile format is the default output format.
%
%   Supports also DICOM image inputs. In DICOM case, you should input a
%   single DICOM image (full path, if the image is not in MATLAB path). ALL
%   images in that folder are then read in the order of the title so you
%   should make sure the images are in the correct numerical/alphabetical
%   order. This feature requires image processing toolbox. For DICOM images
%   you do not need to input pixdim or pixsize. Does not support dynamic
%   frames. DICOM images should have either .dcm or .ima file types.
%
%   Alternatively BMP, PNG or TIFF images can be used. The requirements and
%   usage are the same as in DICOM case, but image processing toolbox is
%   NOT required and pixsize is required. Does not support dynamic frames.
%
% Examples:
%   A = Voxelized_phantom_handle(input_name, output_name, pixdim, pixsize)
%   Voxelized_phantom_handle(input_name, output_name, pixdim, pixsize,
%   first_slice, end_slice, time_frames, lesion_input_file, output_type,
%   windows_format)
%   Voxelized_phantom_handle('phantom_1.bin', 'phantom', [256 256], [0.01
%   0.01 0.01], 'single', 150, 450, 5, 'phantom_lesion_1.bin')
%   Voxelized_phantom_handle('phantom_1.bin', 'phantom', [256 256], [0.01
%   0.01 0.01], [], [], [], [], 'phantom_lesion_1.bin')
%   Voxelized_phantom_handle('image_1.bmp', 'phantom_bmp', [], [0.01
%   0.01 0.01], [], [], [], [], [], 'metaimage')
%   Voxelized_phantom_handle('image_001.dcm', 'phantom_dicom')
%
% Inputs:
%   input_name = Full name of the input data file, e.g. 'phantom_1.bin'.
%   Include full path if the file is not in MATLAB/Octave path. Input file
%   is by default assumed to be in 32-bit float format (images can be in
%   any format supported by the image format itself), but non-float binary
%   data can be used if data_type variable is set accordingly.
%
%   output_name = Name of the output file. The image and header will have
%   the same name, but different file ending (.i33 for image and .h33 for
%   header file when using Interfile and .mhd and .raw when using
%   MetaImage). Use the header file in GATE. 
%
%   pixdim = The image width and height (pixels), e.g. [256 256] if you
%   have 256x256 images. In case of DICOM images, it is taken from the info
%   header and in case of bitmap images, it is taken from the image size.
%   All input images need to be of the same size.
%
%   pixsize = The size (in cm) of each voxel in each ([x y z]) dimension.
%
%   data_type = (optional) The data type of the input data. Default is
%   assumed to 32-bit float. For signed integers use 'int8', 'int16', etc.
%   and for unsigned use 'uint8', 'uint16', etc. For doubles use 'float64'
%   or 'double'.
%
%   first_slice = (optional) Crop the phantom starting from this slice. No
%   further cropping will be performed, even if the specified slice is
%   empty. Can be omitted.
%
%   end_slice = (optional) Crop the phantom ending to this slice. No
%   further cropping will be performed, even if the specified slice is
%   empty. Can be omitted.
%
%   time_frames = (optional) The number of time frames (files) in a dynamic
%   case. Each different time frame/step is assumed to be in a different
%   file with _1, _2, _3, etc. suffix. Can be omitted.
%
%   lesion_input_file = (optional) Full name of the input file containing
%   the lesion data that will be added to the input_name file. Works like
%   input_name and has to have exactly the same properties as the main
%   input file (same sizes, same number of time frames, same file type).
%
%   output_type = (optional) A text string that specifies whether Interfile
%   (default) or MetaImage is used. For Interfile use 'interfile' and for
%   MetaImage use 'metaimage'.
%
%   windows_format = Text file format. Default is Unix style (LF), but
%   setting this value to true uses Windows file format (CR LF). Omitting
%   will use Unix style.
%
% Output:
%   A = The same image (matrix) that was stored in the output file.
%
%
% See also saveInterfile, saveMetaImage, Voxelized_source_handle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Samuli Summala, Ville-Veikko Wettenhovi
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

if nargin >= 5 && ~isempty(varargin) && ~isempty(varargin{1})
    data_type = varargin{1};
else
    data_type = 'single';
end
if nargin >= 6
    first_slice = varargin{2};
    end_slice = varargin{3};
end
if nargin >= 8 && ~isempty(varargin{4})
    tt = varargin{4};
else
    tt = 1;
end
ind = strfind(input_name,'.');
f_type = input_name(ind(end) + 1:end);
if nargin >= 9 && ~isempty(varargin{5})
    lesion_input_file = varargin{5};
    if strcmpi(f_type, 'dcm') || strcmpi(f_type, 'ima')
        f_path = fileparts(which(lesion_input_file));
        if isempty(f_path)
            f_path = fileparts(input_name);
        end
        S = dir(fullfile(f_path,'*.dcm'));
        if isempty(S)
            S = dir(fullfile(f_path,'*.ima'));
        end
        [~, reindex] = sort( str2double( regexp( {S.name}, '\d+', 'match', 'once' )));
        S = S(reindex);
        F = fullfile(f_path,S(1).name);
        X = dicomread(F);
        A = zeros(size(X),numel(S));
        A(:,:,1) = single(X);
        for k = 2:numel(S)
            F = fullfile(f_path,S(k).name);
            X = dicomread(F);
            A(:,:,k) = single(X);
        end
    elseif strcmpi(f_type, 'bmp') || strcmp(f_type, 'png') || strcmp(f_type, 'tiff') || strcmp(f_type, 'tif')
        f_path = fileparts(which(lesion_input_file));
        if isempty(f_path)
            f_path = fileparts(input_name);
        end
        S = dir(fullfile(f_path,['*.' f_type]));
        [~, reindex] = sort( str2double( regexp( {S.name}, '\d+', 'match', 'once' )));
        S = S(reindex);
        F = fullfile(f_path,S(1).name);
        X = dicomread(F);
        info = dicominfo(F);
        pixdim = size(X);
        pixsize = [info.PixelSpacing(1)/10 info.PixelSpacing(2)/10, info.SliceThickness/10];
        A = zeros([pixdim,numel(S)]);
        A(:,:,1) = single(X);
        for k = 2:numel(S)
            F = fullfile(f_path,S(k).name);
            X = dicomread(F);
            A(:,:,k) = single(X);
        end
    else
        fid = fopen(lesion_input_file);
        B = fread(fid,inf,'single=>single',0,'l');
        fclose(fid);
    end
end
if nargin >= 10 && ~isempty(varargin{6})
    output_type = varargin{6};
else
    output_type = 'interfile';
end
if nargin >= 11 && ~isempty(varargin{7})
    windows_format = varargin{7};
else
    windows_format = false;
end
% if nargin >= 9
%     pix_width = varargin{5};
%     peti = varargin{6};
%     pix_space = varargin{7};
% end

if strcmpi(f_type, 'dcm') || strcmpi(f_type, 'ima')
    f_path = fileparts(which(input_name));
    if isempty(f_path)
        f_path = fileparts(input_name);
    end
    S = dir(fullfile(f_path,'*.dcm'));
    if isempty(S)
        S = dir(fullfile(f_path,'*.ima'));
    end
    % source:
    % https://www.mathworks.com/matlabcentral/answers/360531-how-do-i-sort-filenames-containing-text-and-numbers-in-numerical-order-in-matlab
    [~, reindex] = sort( str2double( regexp( {S.name}, '\d+', 'match', 'once' )));
    S = S(reindex);
    %
    F = fullfile(f_path,S(1).name);
    X = dicomread(F);
    info = dicominfo(F);
    pixdim = size(X);
    pixsize = [info.PixelSpacing(1)/10 info.PixelSpacing(2)/10, info.SliceThickness/10];
    A = zeros([pixdim,numel(S)]);
    A(:,:,1) = single(X);
    for k = 2:numel(S)
        F = fullfile(f_path,S(k).name);
        X = dicomread(F);
        A(:,:,k) = single(X);
    end
    tt = 1;
elseif strcmpi(f_type, 'bmp') || strcmp(f_type, 'png') || strcmp(f_type, 'tiff') || strcmp(f_type, 'tif')
    f_path = fileparts(which(input_name));
    if isempty(f_path)
        f_path = fileparts(input_name);
    end
    S = dir(fullfile(f_path,['*.' f_type]));
    [~, reindex] = sort( str2double( regexp( {S.name}, '\d+', 'match', 'once' )));
    S = S(reindex);
    F = fullfile(f_path,S(1).name);
    X = imread(F);
    pixdim = size(X);A = zeros([pixdim,numel(S)]);
    A(:,:,1) = single(X);
    for k = 2:numel(S)
        F = fullfile(f_path,S(k).name);
        X = imread(F);
        A(:,:,k) = single(X);
    end
    tt = 1;
else
    fid = fopen(input_name);
    A = fread(fid,inf,[data_type '=>single'],0,'l');
    fclose(fid);
end
% Add the lesion values to the original values
if nargin >= 9 && ~isempty(varargin{5})
    A = A + B;
end
if strcmpi(f_type, 'dcm') || strcmpi(f_type, 'ima') || strcmpi(f_type, 'bmp') || strcmp(f_type, 'png') || strcmp(f_type, 'tiff') || strcmp(f_type, 'tif')
else
    A = rot90(reshape(A,pixdim(1),pixdim(2),length(A)/(pixdim(1)*pixdim(2))),2);
end
if min(A(:)) < 0
    A = A + abs(min(A(:)));
end
kerroin = 1000 / max(A(:));
if nargin >= 6 && ~isempty(first_slice) && ~isempty(end_slice)
    slices = end_slice - first_slice + 1;
    A = uint16(kerroin*A(:,:,first_slice:end_slice));
else
    slices = size(A,3);
    A = uint16(kerroin*A);
end

cols = reshape(any(A),pixdim(2),slices);
rows = reshape(any(permute(A,[2 1 3])),pixdim(1),slices);
row1 = find(any(rows,2),1,'first');
row2 = find(any(rows,2),1,'last');
col1 = find(any(cols,2),1,'first');
col2 = find(any(cols,2),1,'last');
if nargin < 6 || (nargin >= 6 && isempty(first_slice) && isempty(end_slice))
    slice = squeeze(any(A));
    slice1 = find(any(slice,1),1,'first');
    slice2 = find(any(slice,1),1,'last');
else
    slice1 = 1;
    slice2 = slices;
end

if tt == 1
    A = A(row1:row2,col1:col2,slice1:slice2);
else
    % Leave a slight gap in the dynamic case such that the phantom can
    % freely move in each direction without getting cut-off
    gap = round(max(pixdim) / 50);
    row1 = row1 - gap;
    row2 = row2 + gap;
    col1 = col1 - gap;
    col2 = col2 + gap;
    slice1 = slice1 - gap;
    slice2 = slice2 + gap;
    A = A(row1:row2,col1:col2,slice1:slice2);
end

koko = size(A);
prop.FOV_x = (pixsize(1) * 10) * koko(1);
prop.FOV_y = (pixsize(2) * 10) * koko(2);
prop.axial_FOV = (pixsize(3) * 10) * koko(3);

% Bed support
% Currently unimplemented
% if nargin >= 9
%     gap = ceil((375-334)*pix_space(1)/pix_width)+5;
%     A = [zeros(size(A,1),gap,size(A,3)) A];
%     peti = rot90(imresize(peti,pix_space(1)/pix_width),2);
%     peti = [zeros(floor((size(A,1)-size(peti,1))/2),size(peti,2));peti;zeros(ceil((size(A,1)-size(peti,1))/2),size(peti,2))];
%     peti = [peti zeros(size(A,1),size(A,2)-size(peti,2))];
%     peti(peti>0) = 300;
%     A = A+uint16(peti);
%     disp(['Move the phantom: /gate/phantom/placement/setTranslation 0. ' num2str(-gap*pix_width/2) '. 0. mm'])
% end

if tt == 1
    if strcmpi(output_type, 'metaimage')
        saveMetaImage(output_name, A, [], prop, windows_format)
    else
        saveInterfile(output_name, A, [], [], prop, windows_format)
    end
else
    if strcmpi(output_type, 'metaimage')
        saveMetaImage([output_name '_1'], A, [], prop, windows_format)
    else
        saveInterfile([output_name '_1'], A, [], [], prop, windows_format)
    end
    input_name = input_name(1:ind(end) - 1);
    if nargin >= 8
        lesion_input_file = lesion_input_file(1:strfind(lesion_input_file,'.') - 1);
    end
    if any(strfind(input_name,'_1'))
        input_name = input_name(1:strfind(input_name,'_1') - 1);
        if nargin >= 8
            lesion_input_file = lesion_input_file(1:strfind(lesion_input_file,'_1') - 1);
        end
    end
end

for ll=2:tt
    if nargin >= 8
        fid = fopen([lesion_input_file '_' num2str(ll) '.' f_type]);
        B = fread(fid,inf,'single=>single',0,'l');
        fclose(fid);
    end
    fid = fopen([input_name '_' num2str(ll) '.' f_type]);
    A = fread(fid,inf,'single=>single',0,'l');
    fclose(fid);
    % Add the lesion values to the original values
    if nargin >= 8
        A = A + B;
    end
    A = rot90(reshape(A,pixdim(1),pixdim(2),length(A)/(pixdim(1) * pixdim(2))),2);
    
    if nargin >= 5
        A = uint16(10000*A(row1:row2,col1:col2,first_slice:end_slice));
        A = A(:,:,slice1:slice2);
    else
        A = uint16(10000*A(row1:row2,col1:col2,slice1:slice2));
    end
%     if nargin >= 8
%         A = [zeros(size(A,1),gap,size(A,3)) A];
%         A = A+uint16(peti);
%     end
    if strcmpi(output_type, 'metaimage')
        saveMetaImage([output_name '_frame_' num2str(ll)], A, [], prop, windows_format)
    else
        saveInterfile([output_name '_frame_' num2str(ll)], A, [], [], prop, windows_format)
    end
end
end


