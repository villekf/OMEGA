function options = loadSPECTInterfile(options)
%LOADINTERFILE Loads the input file as a MATLAB matrix with the same
%   dimensions and format type as the original.
%
%   Header file must exist for the file load to work and must contain the
%   number format, number of bytes per pixel and the each individual matrix
%   dimensions. Up to five dimensions are supported natively.
%
%   Supported file types are hdr, img, h33 and i33.
%
% Example:
%   loadInterfile(filename)
%
% Input:
%   filename = Name of the image (.img or .i33) or header file (.hdr or
%   .h33)
%
% See also loadMetaImage
fid = fopen(options.fpath);
if fid == -1
    warning('Specified file was not found! Please select SPECT header file')
    [file, fpath] = uigetfile({'*.*'},'Select SPECT header file');
    if isequal(file, 0)
        error('No file was selected')
    end
    fid = fopen([fpath filesep file]);
end
hdr = textscan(fid,'%s','Delimiter','=');
hdr = hdr{1};
fclose(fid);
% machinefmt = 'l';

ind = find(~cellfun(@isempty, strfind(hdr,'!number of projections'))) + 1;
options.nProjections = str2double(hdr{ind});
ind = find(~cellfun(@isempty, strfind(hdr,'start angle'))) + 1;
options.startAngle = str2double(hdr{ind});
ind = find(~cellfun(@isempty, strfind(hdr,'number of detector heads'))) + 1;
options.nHeads = str2double(hdr{ind});
ind = find(~cellfun(@isempty, strfind(hdr,'extent of rotation'))) + 1;
extRot = str2double(hdr{ind});
options.angleIncrement = extRot / options.nProjections;

options.radiusPerProj = zeros(options.nProjections, 1);
uu = 1;
for kk = 1 : options.nProjections
    ind = find(~cellfun(@isempty, strfind(hdr,['Radius Per View [' num2str(uu) ']'])));
    % ind = find(strfind(hdr{kk},['Radius Per View [' num2str(uu) '] :']));
    options.radiusPerProj(uu) = str2double(hdr{ind + 1});
    uu = uu + 1;
end

% output = fread(fid, inf, [type '=>' type], 0, machinefmt);
% try
%     output = reshape(output, n_dim1, n_dim2, n_dim3, n_dim4, n_dim5);
% catch
%     output = reshape(output, n_dim1, n_dim2, []);
% end
% fclose(fid);