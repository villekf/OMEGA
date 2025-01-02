function [projData, varargout] = loadMCGPUData(size,varargin)
%LOADMCGPUDATA Automatically loads MCGPU projection data
%   This function automatically loads MCGPU projection data from a raw file
%   or files. The user will be prompted for the binary data file when ran.
%   By default only the total image is output, but it is possible to
%   optionally output primary, Compton, Rayleigh and multiple scatter as
%   well. The size of a single projection image needs to be input. Binning
%   can be optionally used.
%
%   Examples:
%       projData = loadMCGPUData([nRowsD nColsD])
%       projData = loadMCGPUData([nRowsD nColsD], binning)
%       projData = loadMCGPUData([nRowsD nColsD], binning, saveData,
%       fileName);
%       projData = loadMCGPUData([nRowsD nColsD], binning, saveData,
%       fileName, filePath);
%       [projData, primaries, Compton, Rayleigh, multipleScatter] =
%       loadMCGPUData([nRowsD nColsD])
%   Inputs:
%       size = See below.
%       nRowsD = The number of rows in the (final) projection image
%       nColsD = The number of columns
%       binning = (optional) The binning factor. Note that if you save
%       primaries, Compton, or any other optional output, binning will be
%       applied to them as well!
%       saveData = (Optional) If set to true, will save the projection
%       images in a mat file with the filename specified with fileName.
%       Note that currently the optional outputs (primaries, etc.) are not
%       saved.
%       fileName = (Optional) The name of the mat-file where the projection
%       images are stored. Required if saveData is set to true.
%       filePath = (Optional) Path to the first projection image. If
%       omitted, the user will be prompted for the file.
%
% See also loadProjectionImages
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2024 Ville-Veikko Wettenhovi
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
nRowsD = size(1);
nColsD = size(2);
if nargin >= 2 && ~isempty(varargin) && ~isempty(varargin{1})
    binning = varargin{1};
else
    binning = 1;
end
if nargin >= 3 && ~isempty(varargin) && ~isempty(varargin{2}) && ~isempty(varargin{3})
    saveData = varargin{2};
    fileName = varargin{3};
else
    saveData = false;
end
if nargin >= 4 && ~isempty(varargin) && ~isempty(varargin{4})
    fpath = varargin{4};
    suffix = strfind(fpath,'.');
    suffix = suffix(end);
    suffix = fpath(suffix + 1 : end);
else
    [file, fpath] = uigetfile('*.raw','Select any projection data file');
    if isequal(file, 0)
        error('No file was selected')
    end
    fpath = fullfile(fpath,file);
    suffix = strfind(file,'.');
    suffix = suffix(end);
    suffix = file(suffix + 1 : end);
end
fid = fopen(fpath);
type = 'single=>single';
A = fread(fid,inf,type);
fclose(fid);
[fpath,file,~] = fileparts(fpath);
%     fnames = dir([fpath '*.' suffix]);
fnames = dir(fpath);
index = ~cellfun('isempty',strfind({fnames.name}, suffix));
fnames = fnames(index);
nFiles = length(fnames);
nProjections = nFiles;
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
if nargout > 1
    varargout{1} = zeros([nRowsD,nColsD,nProjections],class(A));
    if nargout > 2
        varargout{2} = zeros([nRowsD,nColsD,nProjections],class(A));
        if nargout > 3
            varargout{3} = zeros([nRowsD,nColsD,nProjections],class(A));
            if nargout > 4
                varargout{4} = zeros([nRowsD,nColsD,nProjections],class(A));
            end
        end
    end
end
for ll = 1 : nFiles
    if nFiles > 1
        file = allFiles{ll};
        fid = fopen(fullfile(fpath,file));
        A = fread(fid,inf,type);
        fclose(fid);
    end
    proj = A(1 : nRowsD * nColsD);
    proj = reshape(proj, nRowsD * binning, nColsD * binning);
    if nargout > 1
        primary = A(nRowsD * nColsD + 1 : nRowsD * nColsD * 2);
        primary = reshape(primary, nRowsD * binning, nColsD * binning);
        if nargout > 2
            compton = A(nRowsD * nColsD * 2 + 1 : nRowsD * nColsD * 3);
            compton = reshape(compton, nRowsD * binning, nColsD * binning);
            if nargout > 3
                rayleigh = A(nRowsD * nColsD * 3 + 1 : nRowsD * nColsD * 4);
                rayleigh = reshape(rayleigh, nRowsD * binning, nColsD * binning);
                if nargout > 4
                    multiple = A(nRowsD * nColsD * 4 + 1 : end);
                    multiple = reshape(multiple, nRowsD * binning, nColsD * binning);
                end
            end
        end
    end
    if binning > 1
        B = sum(reshape(proj,binning,[]),1,'native');
        B = squeeze(sum(reshape(B,size(A,1)/binning,binning,[]),2,'native'));
        projData(:,:,ll) = B;
        if nargout > 1
            B = sum(reshape(primary,binning,[]),1,'native');
            B = squeeze(sum(reshape(B,size(A,1)/binning,binning,[]),2,'native'));
            varargout{1}(:,:,ll) = B;
            if nargout > 2
                B = sum(reshape(compton,binning,[]),1,'native');
                B = squeeze(sum(reshape(B,size(A,1)/binning,binning,[]),2,'native'));
                varargout{2}(:,:,ll) = B;
                if nargout > 3
                    B = sum(reshape(rayleigh,binning,[]),1,'native');
                    B = squeeze(sum(reshape(B,size(A,1)/binning,binning,[]),2,'native'));
                    varargout{3}(:,:,ll) = B;
                    if nargout > 4
                        B = sum(reshape(multiple,binning,[]),1,'native');
                        B = squeeze(sum(reshape(B,size(A,1)/binning,binning,[]),2,'native'));
                        varargout{4}(:,:,ll) = B;
                    end
                end
            end
        end
    else
        projData(:,:,ll) = proj;
        if nargout > 1
            varargout{1}(:,:,ll) = primary;
            if nargout > 2
                varargout{2}(:,:,ll) = compton;
                if nargout > 3
                    varargout{3}(:,:,ll) = rayleigh;
                    if nargout > 4
                        varargout{4}(:,:,ll) = multiple;
                    end
                end
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