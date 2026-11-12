function saveMetaImage(filename, img, varargin)
%SAVEMETAIMAGE Saves the input matrix as a MetaImage format image
%   This file saves the input 3D or 4D image into MetaImage format. By
%   default the header file is always saved in UNIX text format. You can
%   (optionally) change this into Windows format, regardless of your
%   current OS. This code is based on the code from MathWorks file exchange
%   by Alberto Gomez: 
%   https://se.mathworks.com/matlabcentral/fileexchange/41594-medical-image-processing-toolbox
%
% Examples:
%   saveMetaImage(filename, img)
%   saveMetaImage('image_name', img)
%   saveMetaImage(filename, img, type)
%   saveMetaImage(filename, img, [], options)
%   saveMetaImage(filename, img, [], pixelSize)
%   saveMetaImage(filename, img, type, options, windows_format)
%   saveMetaImage(filename, img, type, pixelSize, windows_format)
%   saveMetaImage(filename, img, type, pixelSize, windows_format, endian)
%
% Input:
%   filename = Name of the image and header files (without file type).
%   Saves the header file as [filename '.mhd'] and image file as [filename
%   '.raw']. Needs to be char input.
%
%   img = The 3D or 4D image. Note that vector format images are not
%   supported.
%
%   type = Data type, e.g. 'single', 'int8', 'uint32', etc. Default is
%   the same type as the input image. If omitted, will use the default
%   value. If the data type is different than the input image type, the
%   image will be casted to the specified new type. For example, if you
%   input int32 image and specify type as int16, the image will be
%   converted to int16 before it is saved.
%
%   options/image_properties = Either the options struct created by the
%   main files or the image_properties struct saved after reconstruction
%   (optional). Necessary if any of the optional values are to be saved
%   (voxel size).
%
%   pixelSize = If saving data where the above struct is not available, this
%   can be used to input the pixel sizes (in mm) for x, y and (possibly) z.
%   This should be input as a vector with the same number of dimensions as
%   the image. Use either the struct or the pixel size, not both. This
%   variable is optional, if omitted values of 1 mm are used by default.
%
%   windows_format = Text file format. Default is UNIX style (LF), but
%   setting this value to true uses Windows file format (CR LF). Omitting
%   will use UNIX style.
%
%   endian = Machine format, i.e. the order of writing bytes. Default is
%   native (usually little endian), but specifying 'b' forces big endian
%   while 'l' forces little endian. 
%
% See also saveImage, saveInterfile

if nargin > 2 && ~isempty(varargin{1})
    type = varargin{1};
else
    type = class(img);
end
if nargin > 3 && ~isempty(varargin{2})
    prop = varargin{2};
else
    prop = [];
end
if nargin > 4 && ~isempty(varargin{3})
    windows_format = varargin{3};
else
    windows_format = false;
end
if nargin > 5 && ~isempty(varargin{4})
    endian = varargin{4};
else
    endian = 'n';
end
img = cast(img,type);
koko = size(img);
if length(koko) == 2
    koko = [koko, 1];
end
if ~isstruct(prop)
    if length(prop) == 2
        prop = [prop, 1];
    end
end
fid = fopen([filename '.raw'],'w');
fwrite(fid, img(:), type, 0, endian);
fclose(fid);
hdrFile = [filename '.mhd'];

writeMetaImageHeader(hdrFile, koko, type, windows_format, endian, prop);
end

function writeMetaImageHeader(hdrFile, koko, type, windows_format, endian, prop)
% The data file is referred to by its name only (no path)
[~, dataName] = fileparts(hdrFile);
dataFile = [dataName '.raw'];
fid = fopen(hdrFile,'w');
element_types=struct('double','MET_DOUBLE','single','MET_FLOAT','int8','MET_CHAR','uint8','MET_UCHAR','int16','MET_SHORT','uint16','MET_USHORT','int32','MET_INT','uint32','MET_UINT',...
    'int64','MET_LONG_LONG', 'uint64', 'MET_ULONG_LONG');

% Number of dimensions (the time dimension of 4D images has a spacing of 1)
nd = length(koko);

%write header file information
fprintf(fid,'ObjectType = Image\r\n');
fprintf(fid,'NDims = %d\r\n', nd);
fprintf(fid,'DimSize = %s\r\n', strtrim(sprintf('%d ', koko(1:nd))));
fprintf(fid,'BinaryData = true\r\n');
if endian == 'b'
    fprintf(fid,'BinaryDataByteOrderMSB = true\r\n');
else
    fprintf(fid,'BinaryDataByteOrderMSB = false\r\n');
end
fprintf(fid,'CompressedData = false\r\n');
fprintf(fid,'TransformMatrix = %s\r\n', strtrim(sprintf('%d ', reshape(eye(nd), 1, []))));
fprintf(fid,'CenterOfRotation = %s\r\n', strtrim(sprintf('%d ', zeros(1, nd))));
fprintf(fid,'ElementType = %s\r\n', element_types.(type));
fprintf(fid,'HeaderSize = 0\r\n');
% Voxel size (mm): either from the field of view of the options/image_properties or from the pixel size vector
voxel = [];
if isstruct(prop)
    % Field of view: FOVa_x/FOVa_y/axial_fov (options), FOVx/FOVy/axialFOV (image_properties) or FOV_x/FOV_y/axial_FOV (old image_properties)
    FOV_x = firstValue(getField(prop, {'FOVa_x', 'FOVx', 'FOV_x'}));
    FOV_y = firstValue(getField(prop, {'FOVa_y', 'FOVy', 'FOV_y'}));
    axial_FOV = firstValue(getField(prop, {'axial_fov', 'axialFOV', 'axial_FOV'}));
    if ~isempty(FOV_x) && ~isempty(FOV_y) && ~isempty(axial_FOV) && FOV_x > 0 && FOV_y > 0 && axial_FOV > 0
        voxel = [FOV_x / koko(1), FOV_y / koko(2), axial_FOV / koko(3)];
    else
        % The voxel sizes saved in image_properties
        voxelX = firstValue(getField(prop, {'voxelSizeX'}));
        voxelY = firstValue(getField(prop, {'voxelSizeY'}));
        voxelZ = firstValue(getField(prop, {'voxelSizeZ'}));
        if ~isempty(voxelX) && ~isempty(voxelY) && ~isempty(voxelZ)
            voxel = [voxelX, voxelY, voxelZ];
        end
    end
elseif isnumeric(prop) && numel(prop) >= 3
    voxel = double(prop(1:3));
end
if isempty(voxel)
    voxel = [1 1 1];
end
voxel = [voxel, ones(1, nd - numel(voxel))];
fprintf(fid,'ElementSpacing = %s\r\n', strtrim(sprintf('%.15g ', voxel(1:nd))));
fprintf(fid,'ElementNumberOfChannels = 1\r\n');
if endian == 'b'
    fprintf(fid,'ElementByteOrderMSB = True\r\n');
else
    fprintf(fid,'ElementByteOrderMSB = False\r\n');
end

fprintf(fid,'ElementDataFile = %s\r\n', dataFile);
fclose(fid);

% The header was written with CR LF line endings. Convert to the requested format (LF by default)
fid = fopen(hdrFile,'r');
fcontent = fread(fid,'uint8');
fclose(fid);
fcontent(fcontent == 13) = [];
if windows_format
    fcontent = uint8(strrep(char(fcontent(:)'), newline, [char(13) newline]))';
end
fid = fopen(hdrFile,'w');
fwrite(fid,fcontent,'uint8');
fclose(fid);

end

function val = getField(prop, names)
%GETFIELD Returns the value of the first of the fields in names that exists in the struct prop, [] if none exist
val = [];
if ~isstruct(prop)
    return
end
for kk = 1 : numel(names)
    if isfield(prop, names{kk})
        val = prop.(names{kk});
        return
    end
end
end

function val = firstValue(x)
%FIRSTVALUE Returns the first element of x as a double (empty if x is empty)
if isempty(x) || (~isnumeric(x) && ~islogical(x))
    val = [];
else
    val = double(x(1));
end
end