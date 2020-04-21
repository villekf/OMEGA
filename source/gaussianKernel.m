function gaussK = gaussianKernel(x, y, z, sigma_x, sigma_y, sigma_z)
%GAUSSIANKERNEL Implements a 3D Gaussian (blurring) kernel
%   Creates a Gaussian kernel from the input values. E.g. if x = -3 : 3,
%   then the x-direction will contain 7 pixels/voxels. Each vector needs to
%   be oriented according to their dimension, e.g. z needs to be 1x1xNz
%   vector and y 1xNyx1 vector. Sigma is the standard deviation. Output is
%   normalized.
%
% Example:
%   gaussK = gaussianKernel(x, y, z, sigma);

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

gaussK = exp(-bsxfun(@plus, bsxfun(@plus, x.^2 / (2*sigma_x^2), y.^2 / (2*sigma_y^2)), z.^2 / (2*sigma_z^2)));

gaussK = gaussK ./  sum(gaussK(:));