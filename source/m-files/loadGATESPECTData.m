function [varargout] = loadGATESPECTData(varargin)
%LOADGATESPECTDATA Loads GATE SPECT projection images
%   Loads, up to 1000, GATE simulated SPECT projection images automatically.
%   Examples:
%       projData = loadGATESPECTData()
%       options = loadGATESPECTData(options)
%       [options, projData] = loadGATESPECTData(options)
%   Inputs:
%       options = The options-struct from a main-file (optional)
%   Outputs:
%       projData = The projection images as a matrix. Each slice (third
%       dimension) corresponds to a specific projection image. If multiple
%       heads are used, the data is loaded as GATE stores them, i.e. first
%       is head 1, then head 2, etc. The output is of type 16-bit unsigned
%       integer (uint16).
%       options = Optionally saves the size of the projection images, the
%       number of projections and the projection images in the options
%       struct.
%
% See also loadInterfile
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2022 Ville-Veikko Wettenhovi
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
if isempty(varargin) || isempty(varargin{1}) || isempty(varargin{1}.fpath)
    [file, fpath] = uigetfile('*.sin','Select the first SPECT projection data file');
    if isequal(file, 0)
        error('No file was selected')
    end
else
    if exist(varargin{1}.fpath,'file') == 0
        warning('Specified file was not found! Please select first projection image')
        [file, fpath] = uigetfile('*.sin','Select the first SPECT projection data file');
        if isequal(file, 0)
            error('No file was selected')
        end
    else
        [fpath,file,~] = fileparts(varargin{1}.fpath);
    end
end
nimi = [fpath file];
nimi = [nimi(1:end-4), '.hdr'];
proj = cell(1000,1);
for kk = 1 : 1000
    try
        proj{kk} = loadInterfile(nimi);
    catch
        break
    end
    if kk < 10
        nimi = [nimi(1:end-5), num2str(kk + 1) '.hdr'];
    else
        nimi = [nimi(1:end-6), num2str(kk + 1) '.hdr'];
    end
end
loppu = kk - 1;
proj(cell2mat(cellfun(@isempty, proj, 'UniformOutput', false))) = [];
koko = cell2mat(cellfun(@size, proj, 'UniformOutput', false));
koko3 = max(koko(:,3));
A = zeros(koko(1,1), koko(1,2), koko3, class(proj{1}));
for kk = 1 : loppu
    cKoko = koko(kk,3);
    A(:,:,1:cKoko) = A(:,:,1:cKoko) + proj{kk};
end
if nargin > 0
    varargout{1} = varargin{1};
    varargout{1}.nRowsD = size(A,1);
    varargout{1}.nColsD = size(A,2);
%     varargout{1}.nProjections = size(A,3);
end
if nargout >= 2
    varargout{2} = A;
    varargout{1}.SinM = A;
elseif nargin == 0
    varargout{1} = A;
else
    varargout{1}.SinM = A;
end
end