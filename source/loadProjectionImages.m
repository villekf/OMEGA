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
%       projData = loadProjectionImages(nProjections, binning, saveData,
%       fileName);
%   Inputs:
%       nProjections = The total number of images (projections) that are
%       loaded. 
%       binning = (Optional) The binning factor.
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
if nargin >= 4 && ~isempty(varargin) && ~isempty(varargin{2}) && ~isempty(varargin{3})
    saveData = varargin{2};
    fileName = varargin{3};
else
    saveData = false;
end
[file, fpath] = uigetfile({'*.tif;*.tiff;*.bmp'},'Select first projection image');
if isequal(file, 0)
    error('No file was selected')
end
t = strfind(file,'_');
if isempty(t)
    t = strfind(file,'0');
    t = t(1) - 1;
else
    t = t(end);
end
suffix = strfind(file,'.');
suffix = suffix(end);
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
        B = sum(reshape(A,binning,[]),1,'native');
        if max(B(:)) == 2^16-1
            B = mean(reshape(A,binning,[]),1);
            B = uint16(squeeze(mean(reshape(B,size(A,1)/binning,binning,[]),2)));
        else
            B = squeeze(sum(reshape(B,size(A,1)/binning,binning,[]),2,'native'));
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