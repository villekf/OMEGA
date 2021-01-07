function vec = deblur(vec, varargin)
%DEBLUR Computes the deblur phase for the PSF reconstruction for the input
%image/vector. Performs symmetric padding.
%
% Examples:
%   vec = deblur(vec, options, gaussK, Nx, Ny, Nz);
%   vec = deblur(vec, g_dim_x, g_dim_y, g_dim_z, gaussK, Nx, Ny, Nz, deblur_iterations);
%
% INPUTS:
%   vec = The image estimes
%   g_dim_x = Size of the PSF kernel in x-direction
%   g_dim_y = Size of the PSF kernel in y-direction
%   g_dim_z = Size of the PSF kernel in z-direction
%   options = Size of PSF kernel (options.g_dim_x, options.g_dim_y and
%   options.g_dim_z), number of deblurring steps
%   (options.deblur_iterations)
%   gaussK = The Gaussian kernel
%   Nx/y/z = Image size in x/y/z direction
%   deblur_iterations = Number of times the deblurring step is performed
%
% OUTPUTS:
%   vec = The deblurred estimates for each iteration
%
% See also gaussianKernel, PSFKernel, computeConvolution

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

if nargin >= 9
    g_dim_x = varargin{1};
    g_dim_y = varargin{2};
    g_dim_z = varargin{3};
    gaussK = varargin{4};
    Nx = varargin{5};
    Ny = varargin{6};
    Nz = varargin{7};
    deblur_iterations = varargin{8};
    epps = 1e-8;
else
    g_dim_x = varargin{1}.g_dim_x;
    g_dim_y = varargin{1}.g_dim_y;
    g_dim_z = varargin{1}.g_dim_z;
    gaussK = varargin{2};
    Nx = varargin{3};
    Ny = varargin{4};
    Nz = varargin{5};
    epps = varargin{1}.epps;
    deblur_iterations = varargin{1}.deblur_iterations;
end
if size(vec,2) == 1
    jelppi = reshape(vec, Nx, Ny, Nz);
    apu = reshape(vec, Nx, Ny, Nz);
else
    jelppi = vec;
    apu = vec;
end
if size(gaussK,2) == 1
    gaussK = reshape(gaussK, g_dim_x*2 + 1, g_dim_y*2 + 1, g_dim_z*2 + 1);
end
apu = padding(apu, [g_dim_x g_dim_y g_dim_z]);
for kk = 1 : deblur_iterations
    apu2 = convn(padding(jelppi, [g_dim_x g_dim_y g_dim_z]), gaussK, 'valid') + epps;
    apu2 = padding(apu2, [g_dim_x g_dim_y g_dim_z]);
    jelppi = jelppi .* convn(apu ./ apu2, gaussK, 'valid');
end
vec = jelppi(:);