function A = Voxelized_source_handle(input_name, output_name, pixdim, pixsize, varargin)
% VOXELIZED_SOURCE_HANDLE This code can be used to crop a voxelized
% source and automatically convert it to uint16 format in Interfile or
% MetaImage file format.
%
%   A voxelized source can be e.g. XCAT or MOBY source. This function
%   automatically crops the source images such that there are no empty
%   rows/columns/slices. Supports both static and dynamic data, though with
%   dynamic data there is slight air gaps in the image to make sure there
%   is no cut-off due to motion. The source is saved in Interfile or
%   MetaImage format with automatically created header file that can be
%   input into GATE/GGEMS. This is practically identical with the
%   Voxelized_phantom_handle, but no scaling will be done and lesions will
%   completely replace the voxels they occupy, i.e. the lesion values
%   higher than zero will replace the corresponding voxel values in the
%   original activity image.
%
%   Like the phantom version, this function also supports both DICOM and
%   BMP/PNG/TIFF formats with the same requirements and restrictions.
%
% Examples:
%   Voxelized_source_handle(input_name, output_name, pixdim, pixsize)
%   Voxelized_source_handle(input_name, output_name, pixdim, pixsize,
%   first_slice, end_slice, row_numbers, col_numbers, time_frames,
%   lesion_input_file, output_type)
%   Voxelized_source_handle('act_1.bin', 'act', [256 256], [0.01
%   0.01 0.01], 150, 450, [5,250], [5,250], 5, 'phantom_lesion_1.bin')
%   Voxelized_source_handle('act_1.bin', 'act', [256 256], [0.01
%   0.01 0.01], [], [], [], [], [], 'act_lesion_1.bin')
%   Voxelized_source_handle('act_1.bin', 'act', [256 256], [0.01
%   0.01 0.01], [], [], [], [], [], [], 'metaimage')
%   Voxelized_source_handle(image, 'act', [256 256], [0.01
%   0.01 0.01], [], [], [], [], [], [], 'metaimage')
%
% Inputs:
%   input_name = Full name of the input data file, e.g. 'act_1.bin'.
%   Include full path if the file is not in MATLAB/Octave path. Binary
%   input file is always assumed to be in 32-bit float format while images
%   can be in any format supported by the image format itself.
%   Alternatively, you can input a 3D grayscale volume instead. In such a
%   case, simply input the volume variable name.
%
%   output_name = Name of the output file. The image and header will have
%   the same name, but different file ending (.i33 for image and .h33 for
%   header file when using Interfile and .mhd and .raw when using
%   MetaImage). Use the header file in GATE.
%
%   pixdim = The image width and height (pixels), e.g. [256 256].
%
%   pixsize = The size (in cm) of each voxel in each ([x y z]) dimension.
%
%   first_slice = (optional) Crop the source starting from this slice. No
%   further cropping will be performed, even if the specified slice is
%   empty. Can be omitted.
%
%   end_slice = (optional) Crop the source ending to this slice. No
%   further cropping will be performed, even if the specified slice is
%   empty. Can be omitted.
%
%   row_numbers = (optional) Crop the source starting from the specified
%   row value and ending to the specified row value. Input should be a 2D
%   vector, for example [x y]. Note that rows in this case mean the rows in
%   MATLAB/Octave (i.e. 1st dimension).
%
%   col_numbers = (optional) Same as above, but for columns.
%
%   time_frames = (optional) The number of time frames (files) in a dynamic
%   case. Each different time frame/step is assumed to be in a different
%   file with _1, _2, _3, etc. suffix.
%
%   lesion_input_file = (optional) Full name of the input file containing
%   the lesion data that will replace the corresponding non-zero voxels in
%   input_name file. Works like input_name and has to have exactly the same
%   properties as the main input file (same sizes, same number of time
%   frames, same file type).
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
% See also saveInterfile, saveMetaImage, Voxelized_phantom_handle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2018-2024 Samuli Summala, Ville-Veikko Wettenhovi
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

if nargin >= 5
    first_slice = varargin{1};
    end_slice = varargin{2};
end
if nargin >= 7 && ~isempty(varargin{3})
    row1 = varargin{3}(1);
    row2 = varargin{3}(2);
else
    row1 = 0;
    row2 = 0;
end
if nargin >= 8 && ~isempty(varargin{4})
    col1 = varargin{4}(1);
    col2 = varargin{4}(2);
else
    col1 = 0;
    col2 = 0;
end
if nargin >= 9
    tt = varargin{5};
    if isempty(tt)
        tt = 1;
    end
else
    tt = 1;
end
if ~isnumeric(input_name)
    ind = strfind(input_name,'.');
    f_type = input_name(ind(end) + 1:end);
    if nargin >= 10 && ~isempty(varargin{6})
        lesion_input_file = varargin{6};
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
            B = zeros(size(X),numel(S));
            B(:,:,1) = single(X);
            for k = 2:numel(S)
                F = fullfile(f_path,S(k).name);
                X = dicomread(F);
                B(:,:,k) = single(X);
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
            X = imread(F);
            pixdim = size(X);
            B = zeros([pixdim,numel(S)]);
            B(:,:,1) = single(X);
            for k = 2:numel(S)
                F = fullfile(f_path,S(k).name);
                X = imread(F);
                B(:,:,k) = single(X);
            end
        else
            fid = fopen(lesion_input_file);
            B = fread(fid,inf,'single=>single',0,'l');
            fclose(fid);
        end
    end
end
if nargin >= 11 && ~isempty(varargin{7})
    output_type = varargin{7};
else
    output_type = 'interfile';
end
if nargin >= 12 && ~isempty(varargin{8})
    windows_format = varargin{8};
else
    windows_format = false;
end
if nargin >= 13 && ~isempty(varargin{9})
    crop = varargin{9};
else
    crop = true;
end
% if nargin >= 9
%     pix_width = varargin{5};
%     peti = varargin{6};
%     pix_space = varargin{7};
% end

if ~isnumeric(input_name)
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
    else
        fid = fopen(input_name);
        A = fread(fid,inf,'single=>single',0,'l');
        fclose(fid);
    end
else
    A = input_name(:,:,:,1);
end

% Replace the original activity values with the lesion values
if nargin >= 10 && ~isempty(varargin{6})
    ind = B > 0;
    A(ind) = B(ind);
    %     A = A + B;
end
if isnumeric(input_name)
elseif (strcmpi(f_type, 'dcm') || strcmpi(f_type, 'ima') || strcmpi(f_type, 'bmp') || strcmp(f_type, 'png') || strcmp(f_type, 'tiff') || strcmp(f_type, 'tif'))
else
    A = rot90(reshape(A,pixdim(1),pixdim(2),length(A)/(pixdim(1)*pixdim(2))),2);
end
if min(A(:)) < 0
    A = A + abs(min(A(:)));
end
kerroin = 1;
if nargin >= 5 && ~isempty(first_slice) && ~isempty(end_slice)
    slices = end_slice - first_slice + 1;
    A = uint16(kerroin*A(:,:,first_slice:end_slice));
else
    slices = size(A,3);
    A = uint16(kerroin*A);
end

cols = reshape(any(A),pixdim(2),slices);
rows = reshape(any(permute(A,[2 1 3])),pixdim(1),slices);
if row1 == 0
    row1 = find(any(rows,2),1,'first');
    row2 = find(any(rows,2),1,'last');
end
if col1 == 0
    col1 = find(any(cols,2),1,'first');
    col2 = find(any(cols,2),1,'last');
end
if nargin < 5 || (nargin >= 5 && isempty(first_slice) && isempty(end_slice))
    slice = squeeze(any(A));
    slice1 = find(any(slice,1),1,'first');
    slice2 = find(any(slice,1),1,'last');
else
    slice1 = 1;
    slice2 = slices;
end
if ~crop
    row1 = 1;
    row2 = size(A,1);
    col1 = 1;
    col2 = size(A,2);
    slice1 = 1;
    slice2 = slices;
end
if nargin >= 5 && ~isempty(first_slice)
    sliceD1 = first_slice;
else
    sliceD1 = slice1;
end
if nargin >= 6 && ~isempty(end_slice)
    sliceD2 = end_slice;
else
    sliceD2 = slice2;
end

if tt == 1
    A = A(row1:row2,col1:col2,slice1:slice2);
    disp(['Figure cropped to [' num2str(row1) ':' num2str(row2) ',' num2str(col1) ':' num2str(col2) ',' num2str(sliceD1) ':' num2str(sliceD2) ']'])
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
    sliceD1 = sliceD1 - gap;
    sliceD2 = sliceD2 + gap;
    A = A(row1:row2,col1:col2,slice1:slice2);
    disp(['Figure cropped to [' num2str(row1) ':' num2str(row2) ',' num2str(col1) ':' num2str(col2) ',' num2str(sliceD1) ':' num2str(sliceD2) ']'])
end

koko = size(A);
if numel(koko) == 2
    koko = [koko, 1];
end
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
    if ~isnumeric(input_name)
        input_name = input_name(1:ind(end) - 1);
        if nargin >= 10 && ~isempty(varargin{6})
            lesion_input_file = lesion_input_file(1:strfind(lesion_input_file,'.') - 1);
        end
        if any(strfind(input_name,'_1'))
            input_name = input_name(1:strfind(input_name,'_1') - 1);
            if nargin >= 10 && ~isempty(varargin{6})
                lesion_input_file = lesion_input_file(1:strfind(lesion_input_file,'_1') - 1);
            end
        end
    end
end

for ll=2:tt
    if isnumeric(input_name)
        A = input_name(:,:,:,tt);
    else
        if nargin >= 10 && ~isempty(varargin{6})
            fid = fopen([lesion_input_file '_' num2str(ll) '.' f_type]);
            B = fread(fid,inf,'single=>single',0,'l');
            fclose(fid);
        end
        fid = fopen([input_name '_' num2str(ll) '.' f_type]);
        A = fread(fid,inf,'single=>single',0,'l');
        fclose(fid);
        % Replace the original activity values with the lesion values
        if nargin >= 8
            ind = B > 0;
            A(ind) = B(ind);
        end
        A = rot90(reshape(A,pixdim(1),pixdim(2),length(A)/(pixdim(1) * pixdim(2))),2);
    end
    if nargin >= 5
        A = uint16(A(row1:row2,col1:col2,first_slice:end_slice));
        A = A(:,:,slice1:slice2);
    else
        A = uint16(A(row1:row2,col1:col2,slice1:slice2));
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

x1 = -(pixsize(1) * 10)*size(A,1) / 2;
y1 = -(pixsize(2) * 10)*size(A,2) / 2;
z1 = -(pixsize(3) * 10)*size(A,3) / 2;
fprintf(['Set the source in GATE: \n /gate/source/voxelsource/setPosition ' num2str(x1) ' ' num2str(y1) ' ' num2str(z1) ' mm \n'])
end