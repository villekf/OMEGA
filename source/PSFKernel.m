function [gaussK, options] = PSFKernel(options)
%PSFKERNEL Compute the Gaussian kernel for the PSF
%   

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
if options.use_psf
    dx = options.FOVa_x / options.Nx;
    dy = options.FOVa_y / options.Ny;
    dz = options.axial_fov / options.Nz;
    
    g_pituus_x = ceil(2*(options.FWHM(1) / (2 * sqrt(2 * log(2)))) / dx);
    g_pituus_y = ceil(2*(options.FWHM(2) / (2 * sqrt(2 * log(2)))) / dy);
    g_pituus_z = ceil(2*(options.FWHM(3) / (2 * sqrt(2 * log(2)))) / dz);
    g_x = linspace(-g_pituus_x * dx, g_pituus_x * dx, 2*g_pituus_x + 1)';
    g_y = linspace(-g_pituus_y * dy, g_pituus_y * dy, 2*g_pituus_y + 1);
    g_z = zeros(1,1,g_pituus_z*2+1);
    g_z(1,1,:) = linspace(-g_pituus_z * dz, g_pituus_z * dz, 2*g_pituus_z + 1);
    gaussK = gaussianKernel(g_x, g_y, g_z, options.FWHM(1) / (2 * sqrt(2 * log(2))), options.FWHM(2) / (2 * sqrt(2 * log(2))), options.FWHM(3) / (2 * sqrt(2 * log(2))));
    
    options.g_dim_x = uint32(g_pituus_x);
    options.g_dim_y = uint32(g_pituus_y);
    options.g_dim_z = uint32(g_pituus_z);
    if options.implementation == 2 || options.implementation == 3
        gaussK = single(gaussK(:));
    end
else
    gaussK = single(0);
end
end

