function A = Voxelized_source_handle(input_name, output_name, pixdim, pixsize, varargin)
% VOXELIZED_SOURCE_HANDLE This code can be used to crop a voxelized
% source and automatically convert it to uint16 format in Interfile file
% format.
%
%   A voxelized source can be e.g. XCAT or MOBY source. This function
%   automatically crops the phantoms such that there are no empty
%   rows/columns/slices. Supports both static and dynamic data, though with
%   dynamic data there is slight air gaps in the phantom to make sure there
%   is no cut-off due to motion. The source is saved in Interfile format
%   with automatically created header file that can be input into GATE.
%   This is practically identical with the Voxelized_phantom_handle, but no
%   scaling will be done and lesions will completely replace the voxels
%   they occupy, i.e. the lesion values higher than zero will replace the
%   corresponding voxel values in the original activity image.
%
% Examples:
%   Voxelized_source_handle(input_name, output_name, pixdim, pixsize)
%   Voxelized_source_handle(input_name, output_name, pixdim, pixsize,
%   first_slice, end_slice, time_frames, lesion_input_file)
%   Voxelized_source_handle('act_1.bin', 'act', [256 256], [0.01
%   0.01 0.01], 150, 450, 5, 'phantom_lesion_1.bin')
%   Voxelized_source_handle('act_1.bin', 'act', [256 256], [0.01
%   0.01 0.01], [], [], [], 'act_lesion_1.bin')
%
% Input:
%   input_name = Full name of the input data file, e.g. 'act_1.bin'.
%   Include full path if the file is not in MATLAB/Octave path. Input file
%   is always assumed to be in 32-bit float format.
%
%   output_name = Name of the output file. The image and header will have
%   the same name, but different file ending (.i33 for image and .h33 for
%   header file). Use the header file in GATE.
%
%   pixdim = The image width and height (pixels), e.g. [256 256].
%
%   pixsize = The size (in cm) of each voxel in each ([x y z]) dimension. 
%
%   first_slice = (optional) Crop the phantom starting from this slice. No
%   further cropping will be performed, even if the specified slice is
%   empty.
%
%   end_slice = (optional) Crop the phantom ending to this slice. No
%   further cropping will be performed, even if the specified slice is
%   empty.
%
%   time_frames = (optional) The number of time frames (files) in a dynamic
%   case. Each different time frame/step is assumed to be in a different
%   file with _1, _2, _3, etc. suffix.
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
%   windows_format = Text file format. Default is UNIX style (LF), but
%   setting this value to true uses Windows file format (CR LF). Omitting
%   will use UNIX style.
%   
%
% See also saveInterfile, saveMetaImage, Voxelized_phantom_handle

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

if nargin >= 5
    first_slice = varargin{1};
    end_slice = varargin{2};
end
if nargin >= 7
    tt = varargin{3};
    if isempty(tt)
        tt = 1;
    end
else
    tt = 1;
end
ind = strfind(input_name,'.');
f_type = input_name(ind(end) + 1:end);
if nargin >= 8
    lesion_input_file = varargin{4};
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
if nargin >= 9 && ~isempty(varargin{5})
    output_type = varargin{5};
else
    output_type = 'interfile';
end
if nargin >= 10 && ~isempty(varargin{6})
    windows_format = varargin{6};
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

if nargin >= 8
    ind = B > 0;
    A(ind) = B(ind);
%     A = A + B;
end
if strcmpi(f_type, 'dcm') || strcmpi(f_type, 'ima') || strcmpi(f_type, 'bmp') || strcmp(f_type, 'png') || strcmp(f_type, 'tiff') || strcmp(f_type, 'tif')
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
row1 = find(any(rows,2),1,'first');
row2 = find(any(rows,2),1,'last');
col1 = find(any(cols,2),1,'first');
col2 = find(any(cols,2),1,'last');
if nargin < 5 || (nargin >= 5 && isempty(first_slice) && isempty(end_slice))
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
    if nargin >= 8
        ind = B > 0;
        A(ind) = B(ind);
    end
    A = rot90(reshape(A,pixdim(1),pixdim(2),length(A)/(pixdim(1) * pixdim(2))),2);
    
    if nargin >= 5
        A = uint16(A(row1:row2,col1:col2,first_slice:end_slice));
        A = A(:,:,slice1:slice2);
    else
        A = uint16(A(row1:row2,col1:col2,slice1:slice2));
    end
    if nargin >= 8
        A = [zeros(size(A,1),gap,size(A,3)) A];
        A = A+uint16(peti);
    end
    if strcmpi(output_type, 'metaimage')
        saveMetaImage([output_name '_frame_' num2str(ll)], A, [], prop, windows_format)
    else
        saveInterfile([output_name '_frame_' num2str(ll)], A, [], [], prop, windows_format)
    end
    %     fid = fopen(output_name,'w+');
    %     fwrite(fid,A,'uint16');
    %     fclose(fid);
end

x1 = -(pixsize(1) * 10)*size(A,1) / 2;
y1 = -(pixsize(2) * 10)*size(A,2) / 2;
z1 = -(pixsize(3) * 10)*size(A,3) / 2;
fprintf(['Set the source in: \n /gate/source/voxelsource/setPosition ' num2str(x1) ' ' num2str(y1) ' ' num2str(z1) ' mm \n'])
end