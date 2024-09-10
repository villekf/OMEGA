% function gaussK = gaussianKernel(x, y, z, sigma_x, sigma_y, sigma_z)
function gaussK = gaussianKernel(x, y, varargin)
%GAUSSIANKERNEL Implements a 3D Gaussian (blurring) kernel
%   Creates a Gaussian kernel from the input values. E.g. if x = -3 : 3,
%   then the x-direction will contain 7 pixels/voxels. Each vector needs to
%   be oriented according to their dimension, e.g. z needs to be 1x1xNz
%   vector and y 1xNyx1 vector. Sigma is the standard deviation. Output is
%   normalized.
%
% Example:
%   gaussK = gaussianKernel(x, y, z, sigma_x, sigma_y, sigma_z);

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
if nargin >= 6
    z = varargin{1};
    sigma_x = varargin{2};
    sigma_y = varargin{3};
    sigma_z = varargin{4};
    gaussK = exp(-bsxfun(@plus, bsxfun(@plus, x.^2 / (2*sigma_x^2), y.^2 / (2*sigma_y^2)), z.^2 / (2*sigma_z^2)));
elseif nargin == 4
    sigma_x = varargin{1};
    sigma_y = varargin{2};
    gaussK = exp(-bsxfun(@plus, x.^2 / (2*sigma_x^2), y.^2 / (2*sigma_y^2)));
else
    error('Invalid number of input parameters')
end
if nargin > 6 && varargin{5}
else
    gaussK = gaussK ./  sum(gaussK(:));
end