function vec = computeConvolution(vec, varargin)
%COMPUTECONVOLUTION Computes symmetrically padded 3D convolution for the
%input vector (image)
%   
% Examples:
%   vec = computeConvolution(vec, options, Nx, Ny, Nz, gaussK);
%   vec = computeConvolution(vec, g_dim_x, g_dim_y, g_dim_z, Nx, Ny, Nz, gaussK);
%
% INPUTS:
%   vec = The image estimes
%   g_dim_x = Size of the PSF kernel in x-direction
%   g_dim_y = Size of the PSF kernel in y-direction
%   g_dim_z = Size of the PSF kernel in z-direction
%   options = Size of PSF kernel (options.g_dim_x, options.g_dim_y and
%   options.g_dim_z)
%   Nx/y/z = Image size in x/y/z direction
%   gaussK = The Gaussian kernel
%
% OUTPUTS:
%   vec = The PSF blurred estimate
%
% See also gaussianKernel, PSFKernel, deblur

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

if nargin >= 8
    g_dim_x = varargin{1};
    g_dim_y = varargin{2};
    g_dim_z = varargin{3};
    Nx = varargin{4};
    Ny = varargin{5};
    Nz = varargin{6};
    gaussK = varargin{7};
else
    g_dim_x = varargin{1}.g_dim_x;
    g_dim_y = varargin{1}.g_dim_y;
    g_dim_z = varargin{1}.g_dim_z;
    Nx = varargin{2};
    Ny = varargin{3};
    Nz = varargin{4};
    gaussK = varargin{5};
end
vec = reshape(vec, Nx, Ny, Nz);
gaussK = reshape(gaussK, g_dim_x * 2 + 1, g_dim_y * 2 + 1, g_dim_z * 2 + 1);
vec = padding(vec, [g_dim_x g_dim_y g_dim_z]);
vec = convn(vec, gaussK, 'valid');
vec = vec(:);
end

