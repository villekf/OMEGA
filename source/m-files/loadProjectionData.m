function [projData] = loadProjectionData(type,varargin)
%LOADPROJECTIONDATA Automatically loads projection (binary) data
%   This function automatically loads projection data from a binary file or
%   files. The user will be prompted for the binary data file when ran.
%   The input type (single, double, int32, etc.) is specified with the
%   first input. If size information is input, then only the portion
%   specified by the size variables is taken. It is also possible to input 
%   the number header variables that are automatically removed either from
%   the beginning or from the end. Furthermore, it is also possible to load
%   all binary files with the same file ending (e.g. all .dat files) from
%   the same directory, where each file is assumed to be a single
%   projection (this is also the default behavior). Binning can also be
%   applied, but no other modifications are performed to the data.
%
%   Examples:
%       projData = loadProjectionData(type)
%       projData = loadProjectionData(type, [nRowsD nColsD nProjections])
%       projData = loadProjectionData(type, [nRowsD nColsD nProjections],
%       binning) 
%       projData = loadProjectionData(type, [nRowsD nColsD nProjections],
%       binning, headerBytes)
%       projData = loadProjectionData(type, [nRowsD nColsD nProjections],
%       binning, headerBytes, loadAll)
%       projData = loadProjectionData(type, [nRowsD nColsD nProjections],
%       binning, headerBytes, loadAll, saveData, fileName);
%   Inputs:
%       type = The data type, e.g. 'single', 'double', 'int16', 'uint32',
%       etc.
%       nRowsD = (optional) The number of rows in the (final) projection
%       image
%       nColsD = (optional) The number of columns
%       nProjections = (optional) The number of projections
%       binning = (optional) The binning factor
%       headerBytes = (optional) The number of header numbers removed from
%       the beginning (when headerBytes > 0) or from the end (headerBytes <
%       0)
%       loadAll = (optional) If true, loads all files in the folder with
%       the same file ending. If nProjections is specified, takes only
%       nProjections number of files alphabetically. Default is true.
%       saveData = (Optional) If set to true, will save the projection
%       images in a mat file with the filename specified with fileName.
%       fileName = (Optional) The name of the mat-file where the projection
%       images are stored. Required if saveData is set to true.
%
% See also loadInveonCTData, loadProjectionImages
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2023 Ville-Veikko Wettenhovi
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
if nargin >= 2 && ~isempty(varargin) && ~isempty(varargin{1})
    nRowsD = varargin{1}(1);
    nColsD = varargin{1}(2);
    nProjections = varargin{1}(3);
else
    nRowsD = 0;
    nColsD = 0;
    nProjections = 0;
end
if nargin >= 3 && ~isempty(varargin) && ~isempty(varargin{2})
    binning = varargin{2};
    if binning > 1 && (nRowsD == 0 || nColsD == 0)
        error('You must input the size information to bin the data')
    end
else
    binning = 1;
end
if nargin >= 4 && ~isempty(varargin) && ~isempty(varargin{3})
    headerBytes = varargin{3};
else
    headerBytes = 0;
end
if nargin >= 5  && ~isempty(varargin) && ~isempty(varargin{4})
    loadAll = varargin{4};
else
    loadAll = true;
end
if nargin >= 7 && ~isempty(varargin) && ~isempty(varargin{5}) && ~isempty(varargin{6})
    saveData = varargin{5};
    fileName = varargin{6};
else
    saveData = false;
end
if nargin >= 8 && ~isempty(varargin) && ~isempty(varargin{7})
    fpath = varargin{7};
    file = dir(fpath);
    index = find(cellfun(@(x) x==0, {file.isdir}, 'UniformOutput', 1) > 0, 1, 'first');
    file = file(index).name;
else
    fpath = [];
end
if isempty(fpath)
    [file, fpath] = uigetfile('*.*','Select any projection data file');
    if isequal(file, 0)
        error('No file was selected')
    end
end
fid = fopen(fullfile(fpath,file));
if headerBytes > 0
    fread(fid,headerBytes,'*uint8');
    type = [type '=>' type];
    A = fread(fid,inf,type);
elseif headerBytes < 0
    fileInfo = dir(fullfile(fpath,file));
    fileSize = fileInfo.bytes;
    if all(type == 'uint8') || all(type == 'logical') || all(type == 'int8')
        nBytes = 1;
    elseif all(type == 'uint16') || all(type == 'int16')
        nBytes = 2;
    elseif all(type == 'uint32') || all(type == 'int32') || all(type == 'single')
        nBytes = 4;
    else
        nBytes = 8;
    end
    type = [type '=>' type];
    A = fread(fid,(fileSize + headerBytes) / nBytes,type);
else
    type = [type '=>' type];
    A = fread(fid,inf,type);
end
fclose(fid);
if loadAll
    suffix = strfind(file,'.');
    suffix = suffix(end);
    suffix = file(suffix + 1 : end);
%     fnames = dir([fpath '*.' suffix]);
    fnames = dir(fpath);
    index = ~cellfun('isempty',strfind({fnames.name}, suffix));
    fnames = fnames(index);
    nFiles = length(fnames);
    if nColsD > 0 && nRowsD > 0 && nProjections > 0
        projData = zeros([nRowsD,nColsD,nProjections],class(A));
    else
        projData = cell(nFiles,1);
        % projData = zeros(numel(A)*nFiles/binning^2, 1, class(A));
    end
    allFiles = cell(nFiles,1);
    apuFile = cell(nFiles,1);
    for ll = 1 : nFiles
        allFiles{ll} = fnames(ll).name;
        apu = fnames(ll).name;
        apuFile{ll} = apu(end-(length(suffix) + 4):end);
    end
    [~,I] = sort(apuFile);
    allFiles = allFiles(I);
else
    nFiles = 1;
    if binning > 1
        if nColsD > 0 && nRowsD > 0 && nProjections > 0
            projData = zeros([nRowsD,nColsD,nProjections],class(A));
        % else
        %     projData = zeros(numel(A)/binning^2, 1, class(A));
        end
    end
end
for ll = 1 : nFiles
    if nFiles > 1
        file = allFiles{ll};
        fid = fopen(fullfile(fpath,file));
        A = fread(fid,inf,type);
        fclose(fid);
        if headerBytes > 0
            A = A(1 + headerBytes : end);
        elseif headerBytes < 0
            A = A(1 : end + headerBytes);
        end
    end
    if nColsD > 0 && nRowsD > 0
        if nFiles == 1
            A = reshape(A, nRowsD * binning, nColsD * binning, nProjections);
        else
            A = reshape(A, nRowsD * binning, nColsD * binning);
        end
    end
    if nColsD > 0 && nRowsD > 0
        if nFiles == 1
            if binning > 1
                for kk = 1 : nProjections
                    B = sum(reshape(A(:,:,kk),binning,[]),1,'native');
                    B = squeeze(sum(reshape(B,size(A,1)/binning,binning,[]),2,'native'));
                    projData(:,:,kk) = B;
                end
            else
                projData = A;
            end
        else
            if binning > 1
                B = sum(reshape(A,binning,[]),1,'native');
                B = squeeze(sum(reshape(B,size(A,1)/binning,binning,[]),2,'native'));
                projData(:,:,ll) = B;
            else
                projData(:,:,ll) = A;
            end
        end
    else
        if nFiles == 1
            projData = A;
        else
            projData{ll} = A;
            if ll == nFiles
                projData = cell2mat(projdata);
            end
        end
    end
end
if saveData
    if exist('OCTAVE_VERSION','builtin') == 0
        save(fileName,'projData','-v7.3')
    else
        save(fileName,'projData','-v7')
    end
end