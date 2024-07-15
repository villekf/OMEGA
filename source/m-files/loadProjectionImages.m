function [projData] = loadProjectionImages(nProjections,varargin)
%LOADPROJECTIONIMAGES Automatically loads projection images from TIFF or
%BMP images
%   This function loads the projection data from TIFF or BMP images. When
%   ran, it will prompt for a single projection image. All projection
%   images in the same folder with the same naming pattern will then be
%   loaded. Images are assumed to be numbered as e.g. _0001, _0002, etc.
%   but the number of zeros can be any, only the underline has to be
%   present before the numbers (i.e. _001, _00001, _01, types are all fine.
%   nProjections number of images are loaded (i.e. files from e.g.
%   _0001,...,_nProjections). Optionally, the user can also input a binning
%   factor. The projection images are then binned from both dimensions
%   according to the binning factor (e.g. with binning = 4, an 1024x1536
%   image becomes 256x384). No further modifications are performed.
%   
%   Examples:
%       projData = loadProjectionImages(nProjections);
%       projData = loadProjectionImages(nProjections, binning);
%       projData = loadProjectionImages(nProjections, binning, fpath);
%       projData = loadProjectionImages(nProjections, binning, fpath,
%       saveData, fileName);
%   Inputs:
%       nProjections = The total number of images (projections) that are
%       loaded. 
%       binning = (Optional) The binning factor.
%       fpath = (Optional) The complete path to the first projection image,
%       if omitted the user will be prompted to select the first projection
%       saveData = (Optional) If set to true, will save the projection
%       images in a mat file with the filename specified with fileName.
%       fileName = (Optional) The name of the mat-file where the projection
%       images are stored. Required if saveData is set to true.
%
% See also loadProjectionData, loadInveonCTData
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
    binning = varargin{1};
else
    binning = 1;
end
if nargin >= 4 && ~isempty(varargin) && ~isempty(varargin{3}) && ~isempty(varargin{4})
    saveData = varargin{3};
    fileName = varargin{4};
else
    saveData = false;
end
if nargin < 3 || (isempty(varargin) && isempty(varargin{2}))
    [file, fpath] = uigetfile({'*.tif;*.tiff;*.bmp'},'Select first projection image');
    if isequal(file, 0)
        error('No file was selected')
    end
else
    if exist(varargin{2},'file') == 0
        warning('Specified file was not found! Please select first projection image')
        [file, fpath] = uigetfile({'*.tif;*.tiff;*.bmp'},'Select first projection image');
        if isequal(file, 0)
            error('No file was selected')
        end
    else
        t = strfind(varargin{2},'/');
        if isempty(t)
            t = strfind(varargin{2},'\');
        end
        file = varargin{2}(t(end) + 1 : end);
        fpath = varargin{2}(1:t(end));
    end
end
t = strfind(file,'_');
suffix = strfind(file,'.');
suffix = suffix(end);
if isempty(t) || suffix - t(end) > 6
    for kk = suffix - 1 : -1 : 1
        testi = str2double(file(kk));
        if isnan(testi)
            break
        end
        t = kk;
    end
else
    t = t(end);
end
nN = suffix - (t + 1);
suffix = file(suffix+1:end);
pFile = file(1:t);
ch = [];
for uu = 1 : nN - 1
    ch = [ch,'0'];
end
ch = [ch,'0'];
lFile = [fpath,pFile,ch,'.',suffix];
zeroStart = true;
if exist(lFile,'file') ~= 2
    ch = [ch(1:end-1),'1'];
    lFile = [fpath,pFile,ch,'.',suffix];
    zeroStart = false;
end
A = imread(lFile);
projData = zeros([size(A)/binning,nProjections],class(A));
summa = false;
for kk = 1 : nProjections
    if zeroStart
        ch = num2str(kk - 1);
    else
        ch = num2str(kk);
    end
    for uu = 1 : nN - numel(ch)
        ch = ['0',ch];
    end
    lFile = [fpath,pFile,ch,'.',suffix];
    A = imread(lFile);
    if binning > 1
        if ~summa
            B = sum(reshape(A,binning,[]),1,'native');
        end
        if max(B(:)) == 2^16-1 || summa
            B = sum(reshape(A,binning,[]),1);
            B = uint32(squeeze(sum(reshape(B,size(A,1)/binning,binning,[]),2)));
            summa = true;
            % B = mean(reshape(A,binning,[]),1);
            % B = uint16(squeeze(mean(reshape(B,size(A,1)/binning,binning,[]),2)));
        else
            B = squeeze(sum(reshape(B,size(A,1)/binning,binning,[]),2,'native'));
            if max(B(:)) == 2^16-1
                B = sum(reshape(A,binning,[]),1);
                B = uint32(squeeze(sum(reshape(B,size(A,1)/binning,binning,[]),2)));
                summa = true;
            end
        end
        if isa(B,'uint32')
            if ~isa(projData,'uint32')
                projData = uint32(projData);
            end
        end
        A = B;
    end
    projData(:,:,kk) = A;
end
if saveData
    if exist('OCTAVE_VERSION','builtin') == 0
        save(fileName,'projData','-v7.3')
    else
        save(fileName,'projData','-v7')
    end
end