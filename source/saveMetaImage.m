function saveMetaImage(filename, img, varargin)
%SAVEMETAIMAGE Saves the input matrix as a MetaImage format image
%   This file saves the input 3D or 4D image into MetaImage format (32-bit
%   float by default). This code is based on the code from MathWorks file
%   exchange by Alberto Gomez:
%   https://se.mathworks.com/matlabcentral/fileexchange/41594-medical-image-processing-toolbox
%
% Examples:
%   saveInterfile(filename, img)
%   saveInterfile(filename, img, type)
%   saveInterfile(filename, img, [], options)
%
% Input:
%   filename = Name of the image and header files (without file type).
%   Saves the header file as [filename '.mhd'] and image file as [filename
%   '.raw'].
%
%   img = The 3D or 4D image
%
%   type = Data type, e.g. 'single', 'int8', 'uint32', etc. Default is
%   the same type as the input image.
%
%   options/image_properties = Either the options struct created by the
%   main files or the image_properties struct saved in the cell-matrix
%   (optional). Necessary if any of the optional values are to be saved
%   (voxel size).

if nargin > 3 && ~isempty(varargin{2})
    prop = varargin{2};
else
    prop = [];
end
if nargin > 2 && ~isempty(varargin{1})
    type = varargin{1};
else
    type = class(img);
end
img = cast(img,type);
koko = size(img);
fid = fopen([filename '.raw'],'w');
fwrite(fid, img(:), type);
fclose(fid);
hdrFile = [filename '.mhd'];

if nargin > 3
    writeMetaImageHeader(hdrFile, koko, type, prop);
else
    writeMetaImageHeader(hdrFile, koko, type);
end
end

function writeMetaImageHeader(hdrFile, koko, type, varargin)
if nargin >= 4 && ~isempty(varargin{1})
    prop = varargin{1};
else
    prop = [];
end
fid = fopen(hdrFile,'w');
element_types=struct('double','MET_DOUBLE','int8','MET_CHAR','uint8','MET_UCHAR','int16','MET_SHORT','uint16','MET_USHORT','int32','MET_INT','uint32','MET_UINT',...
    'int64','MET_LONG', 'uint64', 'MET_ULONG');

%write header file information
fprintf(fid,'ObjectType = Image\r\n');
fprintf(fid,['NDims = ' num2str(length(koko)) '\r\n']);
s = 'DimSize =';
for kk = 1 : length(koko)
    s = [s, ' ' num2str(koko(kk))];
end
fprintf(fid,[s '\r\n']);
fprintf(fid,['ElementType = ' element_types.(type) '\r\n']);
fprintf(fid,'HeaderSize = 0\r\n');
if nargin >= 4
    if isfield(prop,'FOV_x')
        fprintf(fid,['ElementSpacing = ' num2str(prop.FOV_x/koko(1)) ' ' num2str(prop.FOV_y/koko(2)) ' ' num2str(prop.axial_FOV/koko(3)) '\r\n']);
    elseif isfield(prop,'FOVa_x')
        fprintf(fid,['ElementSpacing = ' num2str(prop.FOVa_x/koko(1)) ' ' num2str(prop.FOVa_y/koko(2)) ' ' num2str(prop.axial_fov/koko(3)) '\r\n']);
    else
        fprintf(fid,['ElementSpacing = 1 1 1\r\n']);
    end
else
    fprintf(fid,['ElementSpacing = 1 1 1\r\n']);
end
fprintf(fid,'ElementByteOrderMSB = False\r\n');
fprintf(fid,['ElementDataFile = ' hdrFile(1:end-3) 'raw\r\n']);
fclose(fid);

end