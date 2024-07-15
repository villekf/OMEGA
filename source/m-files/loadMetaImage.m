function [output, varargout] = loadMetaImage(filename)
%LOADMETAIMAGE Loads the input file as a MATLAB matrix with the same
%   dimensions and format type as the original.
%
%   Header file must exist for the file load to work and must contain the
%   Number of dimensions, Size of each dimension, Data type and Name of the
%   data file. Header size is used if available. All other data fields are
%   ignored.
%
%   Supported file types are mhd and mha.
%   
%   This code is based on the code from MathWorks file exchange by Alberto
%   Gomez:
%   https://se.mathworks.com/matlabcentral/fileexchange/41594-medical-image-processing-toolbox
%
% Examples:
%   output = loadMetaImage(filename)
%   [output, struct] = loadMetaImage(filename)
%
% Input:
%   filename = Name of the header file (.mhd or .mha)
%
% Outputs:
%   output = The image stored
%   struct = A struct containing many of the possible variables included
%
% See also loadInterfile
element_types = struct('MET_DOUBLE','double','MET_FLOAT', 'single', 'MET_CHAR','int8','MET_UCHAR','uint8','MET_SHORT','int16','MET_USHORT','uint16','MET_INT','int32','MET_UINT','uint32',...
    'MET_LONG', 'int64', 'MET_ULONG','uint64');
n_bytes = struct('double',8,'single', 4, 'int8',1,'uint8',1,'int16',2,'uint16',2,'int32',4,'uint32',4, 'int64',8, 'uint64', 8);
tline = 0;
ll = 1;
M = cell(1,1);
N = cell(1,1);
if strcmp(filename(end-2:end),'raw')
    filename = [filename(1:end-3) 'mhd'];
end
fid = fopen(filename);
if fid < 0
    error(['Could not find ' filename])
end
while sum(tline ~= -1) >= 1 || isempty(tline ~= -1)
    tline = fgetl(fid);
    M{ll} = {lower(tline)};
    N{ll} = {tline};
    ll = ll + 1;
end
fclose(fid);
M = M(1:end-1);
NDims = [];
skip = 0;
type = 'single';
endian = 'l';
transformmatrix = [];

for kk = 1 : length(M)
    apu = cell2mat(M{kk});
    ind = strfind(apu,'=');
    if ~cellfun('isempty',strfind(M{kk},'ndims'))
        NDims = str2double(apu(ind(1) + 1:end));
        if nargout >= 2
            varargout{1}.NDims = NDims;
        end
    elseif ~cellfun('isempty',strfind(M{kk},'dimsize'))
        if isempty(NDims)
            for ii = 1 : length(M)
                apu2 = cell2mat(M{ii});
                ind = strfind(apu2,'=');
                if ~cellfun('isempty',strfind(M{ii},'ndims'))
                    NDims = str2double(apu2(ind(1) + 1:end));
                end
            end
        end
        apu = apu(ind+1:end);
        ind = strfind(apu,' ');
        ind = [ind, length(apu)];
        Dims = zeros(NDims,1);
        for ll = 1 : NDims
            Dims(ll) = str2double(apu(ind(ll):ind(ll + 1)));
        end
        if nargout >= 2
            varargout{1}.Dims = Dims;
        end
    elseif ~cellfun('isempty',strfind(M{kk},'headersize'))
        skip = str2double(apu(ind(1) + 1:end));
        if nargout >= 2
            varargout{1}.headersize = skip;
        end
    elseif ~cellfun('isempty',strfind(M{kk},'elementtype'))
        apu = apu(ind+1:end);
        ind = strfind(apu,' ');
        type = upper(apu(ind + 1:end));
        type = element_types.(type);
        if nargout >= 2
            varargout{1}.element_type = type;
        end
    elseif ~cellfun('isempty',strfind(M{kk},'elementdatafile'))
        apu = cell2mat(N{kk});
        ind = strfind(apu,'=');
        apu = apu(ind+1:end);
        ind = strfind(apu,' ');
        raw_filename = apu(ind + 1:end);
        if any(strfind(raw_filename, '/'))
            raw_filename = raw_filename(strfind(raw_filename, '/') + 1 : end);
        end
        if nargout >= 2
            varargout{1}.raw_filename = raw_filename;
        end
    elseif ~cellfun('isempty',strfind(M{kk},'transformmatrix'))
        apu = apu(ind+1:end);
        ind = strfind(apu,' ');
        transformmatrix = str2num(apu(ind + 1:end));
        if nargout >= 2
            varargout{1}.transformmatrix = transformmatrix;
        end
    elseif ~cellfun('isempty',strfind(M{kk},'elementbyteordermsb')) || ~cellfun('isempty',strfind(M{kk},'binarydatabyteordermsb'))
        apu = apu(ind+1:end);
        ind = strfind(apu,' ');
        if any(strfind(apu(ind + 1:end),'true'))
            endian = 'b';
        end
        if nargout >= 2
            varargout{1}.byteOrderMSB = apu(ind + 1:end);
        end
    elseif ~cellfun('isempty',strfind(M{kk},'elementspacing')) && nargout >= 2
        if isempty(NDims)
            for ii = 1 : length(M)
                apu2 = cell2mat(M{ii});
                ind = strfind(apu2,'=');
                if ~cellfun('isempty',strfind(M{ii},'ndims'))
                    NDims = str2double(apu2(ind(1) + 1:end));
                end
            end
        end
        apu = apu(ind+1:end);
        ind = strfind(apu,' ');
        ind = [ind, length(apu)];
        EleSpacing = zeros(NDims,1);
        for ll = 1 : NDims
            EleSpacing(ll) = str2double(apu(ind(ll):ind(ll + 1)));
        end
        varargout{1}.EleSpacing = EleSpacing;
    elseif ~cellfun('isempty',strfind(M{kk},'offset')) && nargout >= 2
        if isempty(NDims)
            for ii = 1 : length(M)
                apu2 = cell2mat(M{ii});
                ind = strfind(apu2,'=');
                if ~cellfun('isempty',strfind(M{ii},'ndims'))
                    NDims = str2double(apu2(ind(1) + 1:end));
                end
            end
        end
        apu = apu(ind+1:end);
        ind = strfind(apu,' ');
        ind = [ind, length(apu)];
        Offset = zeros(NDims,1);
        for ll = 1 : NDims
            Offset(ll) = str2double(apu(ind(ll):ind(ll + 1)));
        end
        varargout{1}.Offset = Offset;
    end
end

if strcmpi(raw_filename, 'local')
    skip = -1;
    fid = fopen(filename);
else
    D = dir(filename);
    filename = [D.folder '/' raw_filename];
    fid = fopen(filename);
    if fid < 0
        error(['Could not find ' filename])
    end
end
if skip == -1
    D = dir(filename);
    skip = D.bytes - n_bytes.(type) * prod(Dims) ;
    fread(fid, skip, '*uint8');
end
output = fread(fid, prod(Dims), [type '=>' type], endian);
output = reshape(output, Dims');
% if isempty(transformmatrix)
%     transformmatrix = reshape(eye(NDims), 1,NDims*NDims);
% else
%     transformmatrix = reshape(transformmatrix, NDims, NDims);
%     output = reshape(output, NDims, numel(output)/NDims);
%     output = cast(double(transformmatrix) * double(output), type);
%     output = reshape(output, Dims');
% end
fclose(fid);