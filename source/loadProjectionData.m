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
%   projection. Binning can also be applied, but no other modifications are
%   performed to the data.
%
%   Examples:
%       projData = loadProjectionData(type)
%       projData = loadProjectionData(type, [ySize xSize nProjections])
%       projData = loadProjectionData(type, [ySize xSize nProjections],
%       binning) 
%       projData = loadProjectionData(type, [ySize xSize nProjections],
%       binning, headerBytes)
%       projData = loadProjectionData(type, [ySize xSize nProjections],
%       binning, headerBytes, loadAll)
%       projData = loadProjectionData(type, [ySize xSize nProjections],
%       binning, headerBytes, loadAll, saveData, fileName);
%   Inputs:
%       type = The data type, e.g. 'single', 'double', 'int16', 'uint32',
%       etc.
%       ySize = (optional) The number of rows in the projection image
%       xSize = (optional) The number of columns
%       nProjections = (optional) The number of projections
%       binning = (optional) The binning factor
%       headerBytes = (optional) The number of header numbers removed from
%       the beginning (when headerBytes > 0) or from the end (headerBytes <
%       0)
%       loadAll = (optional) If true, loads all files in the folder with
%       the same file ending. If nProjections is specified, takes only
%       nProjections number of files alphabetically.
%       saveData = (Optional) If set to true, will save the projection
%       images in a mat file with the filename specified with fileName.
%       fileName = (Optional) The name of the mat-file where the projection
%       images are stored. Required if saveData is set to true.
%
% See also loadInveonCTData, loadProjectionImages
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021 Ville-Veikko Wettenhovi
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
    ySize = varargin{1}(1);
    xSize = varargin{1}(2);
    nProjections = varargin{1}(3);
else
    ySize = 0;
    xSize = 0;
    nProjections = 0;
end
if nargin >= 3 && ~isempty(varargin) && ~isempty(varargin{2})
    binning = varargin{2};
    if binning > 1 && (ySize == 0 || xSize == 0)
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
    loadAll = false;
end
if nargin >= 7 && ~isempty(varargin) && ~isempty(varargin{5}) && ~isempty(varargin{6})
    saveData = varargin{5};
    fileName = varargin{6};
else
    saveData = false;
end
type = [type '=>' type];
[file, fpath] = uigetfile('*.*','Select any projection data file');
if isequal(file, 0)
    error('No file was selected')
end
fid = fopen(fullfile(fpath,file));
A = fread(fid,inf,type);
fclose(fid);
if headerBytes > 0
    A = A(1 + headerBytes : end);
elseif headerBytes < 0
    A = A(1 : end + headerBytes);
end
if loadAll
    suffix = strfind(file,'.');
    suffix = suffix(end);
    suffix = file(suffix + 1 : end);
    fnames = dir([fpath '*.' suffix]);
    nFiles = length(fnames);
    if xSize > 0 && ySize > 0 && nProjections > 0
        projData = zeros([ySize/binning,xSize/binning,nProjections],class(A));
    else
        projData = zeros(numel(A)*nFiles/binning^2, 1, class(A));
    end
else
    nFiles = 1;
    if binning > 1
        if xSize > 0 && ySize > 0 && nProjections > 0
            projData = zeros([ySize/binning,xSize/binning,nProjections],class(A));
        else
            projData = zeros(numel(A)/binning^2, 1, class(A));
        end
    end
end
for ll = 1 : nFiles
    if ll > 1
        file = fnames(ll).name;
        fid = fopen(fullfile(fpath,file));
        A = fread(fid,inf,type);
        fclose(fid);
        if headerBytes > 0
            A = A(1 + headerBytes : end);
        elseif headerBytes < 0
            A = A(1 : end + headerBytes);
        end
    end
    if xSize > 0 && ySize > 0
        if nFiles == 1
            A = reshape(A, ySize, xSize, nProjections);
        else
            A = reshape(A, ySize, xSize);
        end
    end
    if xSize > 0 && ySize > 0
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
    end
end
if saveData
    if exist('OCTAVE_VERSION','builtin') == 0
        save(fileName,'projData','-v7.3')
    else
        save(fileName,'projData','-v7')
    end
end