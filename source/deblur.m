function vec = deblur(vec, options, ~, ~, gaussK, Nx, Ny, Nz)
%DEBLUR Computes the deblur phase for the PSF reconstruction for the input
%image/vector. Performs symmetric padding.
%
% Example:
%   vec = deblur(vec, options, iter, subsets, gaussK, Nx, Ny, Nz);
%
% INPUTS:
%   vec = The image estimes
%   options = Size of PSF kernel
%   iter = Current iteration number (unused)
%   subset = Total number of subsets (unused)
%   gaussK = The Gaussian kernel
%   Nx/y/z = Image size in x/y/z direction
%
% OUTPUTS:
%   vec = The deblurred estimates for each iteration
%
% See also gaussianKernel, PSFKernel, computeConvolution

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

if size(vec,2) == 1
    jelppi = reshape(vec, Nx, Ny, Nz);
    apu = reshape(vec, Nx, Ny, Nz);
else
    jelppi = vec;
    apu = vec;
end
apu = padding(apu, [options.g_dim_x options.g_dim_y options.g_dim_z]);
for kk = 1 : options.deblur_iterations
    apu2 = convn(padding(jelppi, [options.g_dim_x options.g_dim_y options.g_dim_z]), gaussK, 'valid');
    apu2 = padding(apu2, [options.g_dim_x options.g_dim_y options.g_dim_z]);
    jelppi = jelppi .* convn(apu ./ apu2, gaussK, 'valid');
end
vec = jelppi(:);