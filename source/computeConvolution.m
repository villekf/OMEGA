function vec = computeConvolution(vec, options, Nx, Ny, Nz, gaussK)
%COMPUTECONVOLUTION Computes symmetrically padded 3D convolution for the
%input vector (image)
%   
% Example:
%   vec = computeConvolution(vec, options, Nx, Ny, Nz, gaussK);
%
% INPUTS:
%   vec = The image estimes
%   options = Size of PSF kernel
%   Nx/y/z = Image size in x/y/z direction
%   gaussK = The Gaussian kernel
%
% OUTPUTS:
%   vec = The PSF blurred estimate
%
% See also gaussianKernel, PSFKernel, deblur

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

vec = reshape(vec, Nx, Ny, Nz);
gaussK = reshape(gaussK, options.g_dim_x * 2 + 1, options.g_dim_y * 2 + 1, options.g_dim_z * 2 + 1);
vec = padding(vec, [options.g_dim_x options.g_dim_y options.g_dim_z]);
vec = convn(vec, gaussK, 'valid');
vec = vec(:);
end

