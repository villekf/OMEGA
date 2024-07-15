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

if nargin > 3
    writeMetaImageHeader(hdrFile, koko, type, windows_format, endian, prop);
else
    writeMetaImageHeader(hdrFile, koko, type, windows_format, endian);
end
end

function writeMetaImageHeader(hdrFile, koko, type, windows_format, endian, varargin)
if nargin >= 5 && ~isempty(varargin) && ~isempty(varargin{1})
    prop = varargin{1};
else
    prop = [];
end
fid = fopen(hdrFile,'w');
element_types=struct('double','MET_DOUBLE','single','MET_FLOAT','int8','MET_CHAR','uint8','MET_UCHAR','int16','MET_SHORT','uint16','MET_USHORT','int32','MET_INT','uint32','MET_UINT',...
    'int64','MET_LONG', 'uint64', 'MET_ULONG');

%write header file information
fprintf(fid,'ObjectType = Image\r\n');
fprintf(fid,['NDims = ' num2str(length(koko)) '\r\n']);
s = 'DimSize =';
for kk = 1 : length(koko)
    s = [s, ' ' num2str(koko(kk))];
end
fprintf(fid,[s '\r\n']);
fprintf(fid,'BinaryData = true\n');
if endian == 'b'
    fprintf(fid,'BinaryDataByteOrderMSB = true\n');
else
    fprintf(fid,'BinaryDataByteOrderMSB = false\n');
end
fprintf(fid,'CompressedData = false\n');
fprintf(fid,'TransformMatrix = 1 0 0 0 1 0 0 0 1\n');
fprintf(fid,'CenterOfRotation = 0 0 0\n');
fprintf(fid,['ElementType = ' element_types.(type) '\r\n']);
fprintf(fid,'HeaderSize = 0\r\n');
if nargin >= 4
    if isfield(prop,'FOV_x')
        fprintf(fid,['ElementSpacing = ' num2str(prop.FOV_x/koko(1)) ' ' num2str(prop.FOV_y/koko(2)) ' ' num2str(prop.axial_FOV/koko(3)) '\r\n']);
    elseif isfield(prop,'FOVa_x')
        fprintf(fid,['ElementSpacing = ' num2str(prop.FOVa_x/koko(1)) ' ' num2str(prop.FOVa_y/koko(2)) ' ' num2str(prop.axial_fov/koko(3)) '\r\n']);
    elseif isfield(prop,'FOVx')
        fprintf(fid,['ElementSpacing = ' num2str(prop.FOVx(1)/koko(1)) ' ' num2str(prop.FOVy(1)/koko(2)) ' ' num2str(prop.axialFOV(1)/koko(3)) '\r\n']);
    elseif numel(prop) == numel(koko)
        fprintf(fid,['ElementSpacing = ' num2str(prop(1)) ' ' num2str(prop(2)) ' ' num2str(prop(3)) '\r\n']);
    elseif isfield(prop,'voxelSiceX')
        fprintf(fid,['ElementSpacing = ' num2str(prop.voxelSiceX) ' ' num2str(prop.voxelSiceY) ' ' num2str(prop.voxelSiceZ) '\r\n']);
    else
        fprintf(fid,['ElementSpacing = 1 1 1\r\n']);
    end
else
    fprintf(fid,['ElementSpacing = 1 1 1\r\n']);
end
fprintf(fid,'ElementNumberOfChannels = 1\n');
if endian == 'b'
    fprintf(fid,'ElementByteOrderMSB = True\r\n');
else
    fprintf(fid,'ElementByteOrderMSB = False\r\n');
end

fprintf(fid,['ElementDataFile = ' hdrFile(1:end-3) 'raw\r\n']);
fclose(fid);
if (ispc && ~windows_format) || ((isunix || ismac) && windows_format)
    fid = fopen(hdrFile,'r');
    % This part of code is from: https://www.mathworks.com/matlabcentral/fileexchange/13651-unix2dos
    fcontent = fread(fid,'uint8');
    fclose(fid);
    if ispc && ~windows_format
        CR = char(13);
        fcontent(fcontent==CR) = [];
    elseif (isunix || ismac) && windows_format
        CR = char(13);
        LF = char(10);
        fcontent = strrep(char(row(fcontent)),LF,[CR LF]);
    end
    fid=fopen(hdrFile,'w');
    fwrite(fid,fcontent,'uint8');
    fclose(fid);
end

end