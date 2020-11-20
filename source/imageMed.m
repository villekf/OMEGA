function imageMed(img,varargin)
%IMAGEMED Plots the input 3D image in the specified plane
%   This function can visualize a 3D matrix/image in the desired plane.
%   There are two possible uses, one possibility is inputting the desired
%   plane with the input data. This is achieved by e.g.
%   imageMed(img(2,:,:)) which will then visualize the 2nd and 3rd
%   dimension values on the 2nd row (1st dimension). Second use, is to
%   visualize all three planes of the input image, e.g. imageMed(img, 5)
%   will visualize the 5th planes in all three dimensions. It is also
%   possibly to specify individual planes, e.g. imageMed(img, 5, 4, 2).
%   This function supports pz-cell data, in which case the last iteration
%   is always visualized. Currently dynamic data is not supported. For 4D
%   images, the last slice is taken.
%
%   Examples:
%       imageMed(img(:,5,:))
%       imageMed(img, 5);
%       imageMed(img, 5, 6, 7);
%       imageMed(pz{2}, 5, 6, 7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ville-Veikko Wettenhovi
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
[i, j, k, p, t] = size(img);
if t > 1
    error('Dynamic data is not supported')
end
if nargin == 1
    if i == 1
        apu = permute(img(:,:,:,end), [2 3 1]);
        imagesc(apu)
    elseif j == 1
        apu = permute(img(:,:,:,end), [1 3 2]);
        imagesc(apu)
    elseif k == 1
        imagesc(img)
    end
    axis image
else
    if nargin == 2
        if varargin{1} > i
            error('Second input exceeds the first dimension of the image!')
        end
        apu1 = permute(img(varargin{1},:,:,end), [2 3 1]);
        if varargin{1} > j
            error('Second input exceeds the second dimension of the image!')
        end
        apu2 = permute(img(:,varargin{1},:,end), [1 3 2]);
        if varargin{1} > k
            error('Second input exceeds the third dimension of the image!')
        end
        subplot 131
        imagesc(apu1)
        axis image
        subplot 132
        imagesc(apu2)
        axis image
        subplot 133
        imagesc(img(:,:,varargin{1},end))
        axis image
    elseif nargin == 4
        if varargin{1} > i
            error('Second input exceeds the first dimension of the image!')
        end
        apu1 = permute(img(varargin{1},:,:,end), [2 3 1]);
        if varargin{2} > j
            error('Third input exceeds the second dimension of the image!')
        end
        apu2 = permute(img(:,varargin{2},:,end), [1 3 2]);
        if varargin{3} > k
            error('Fourth input exceeds the third dimension of the image!')
        end
        subplot 131
        imagesc(apu1)
        axis image
        subplot 132
        imagesc(apu2)
        axis image
        subplot 133
        imagesc(img(:,:,varargin{3},end))
        axis image
    else
        error('Unsupported number of input arguments! There has to be either 1, 2 or 4 input arguments.')
    end
end
end