function [gaussK, options, varargout] = PSFKernel(options,varargin)
%PSFKERNEL Compute the Gaussian kernel for the PSF
% Returns the Gaussian PSF kernel and the number of pixels inside the
% specified FWHM. This function can operate either with the options-struct
% or manual inputs. With manual inputs, both 2D and 3D versions are
% available. 2D is used if the number of dimensions of the FWHM is 2.
%
% EXAMPLES:
%   [gaussK, options] = PSFKernel(options)
%   [gaussK, g_dim_x, g_dim_y, g_dim_z] = PSFKernel(Nx, Ny, Nz, FOVa_x FOVa_y, axial_fov, FWHM, implementation)
%   [gaussK, g_dim_x, g_dim_y, g_dim_z] = PSFKernel(FWHM, dx, dy, dz, implementation)
%   [gaussK, g_dim_x, g_dim_y] = PSFKernel(FWHM, dx, dy, implementation)
% INPUTS:
%   Nx = Image size in x-dimension
%   Ny = Image size in y-dimension
%   Nz = Image size in z-dimension
%   FOVa_x = FOV size in x-dimension
%   FOVa_y = FOV size in y-dimension
%   axial_fov = FOV size in z-dimension
%   FWHM = FWHM of the Gaussian kernel in all three dimensions [x y
%   z] when using 3D data, or two dimensions [x y] when using 2D data
%   dx = Pixel size in x-dimension
%   dy = Pixel size in y-dimension
%   dz = Pixel size in z-dimension (not required for 2D data (i.e. can be
%   omitted)
%   options.implementation = Converts the output data to single precision
%   if implementations 2 or 3 are selected (optional)
%
% See also computeConvolution


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
if ~isstruct(options) && nargin >= 7
    Nx = options;
    clear options
    options.use_psf = true;
    options.Nx = Nx;
    options.Ny = varargin{1};
    options.Nz = varargin{2};
    options.FOVa_x = varargin{3};
    options.FOVa_y = varargin{4};
    options.axial_fov = varargin{5};
    options.FWHM = varargin{6};
    if nargin >= 8
        options.implementation = varargin{7};
    end
elseif ~isstruct(options) && nargin >= 3 && ndims(options) >= 2
    FWHM = options;
    clear options
    options.FWHM = FWHM;
    dx = varargin{1};
    dy = varargin{2};
    if numel(options.FWHM) == 3 && options.FWHM(3) > 0
        dz = varargin{3};
        if nargin >= 5
            options.implementation = varargin{4};
        end
    else
        if nargin >= 4
            options.implementation = varargin{3};
        end
    end
end
if ~isfield(options,'use_psf') || options.use_psf
%     if numel(options.FWHM) == 2
%         options.FWHM(3) = 0;
%     end
    if ~exist('dx','var')
        dx = options.FOVa_x / double(options.Nx);
        dy = options.FOVa_y / double(options.Ny);
        if numel(options.FWHM) == 3 && options.FWHM(3) > 0
            dz = options.axial_fov / double(options.Nz);
        end
    end
    
    g_pituus_z = 0;
    g_pituus_x = ceil(2*(options.FWHM(1) / (2 * sqrt(2 * log(2)))) / dx);
    g_pituus_y = ceil(2*(options.FWHM(2) / (2 * sqrt(2 * log(2)))) / dy);
    if numel(options.FWHM) == 3 && options.FWHM(3) > 0
        g_pituus_z = ceil(2*(options.FWHM(3) / (2 * sqrt(2 * log(2)))) / dz);
    end
    g_x = linspace(-g_pituus_x * dx, g_pituus_x * dx, 2*g_pituus_x + 1)';
    g_y = linspace(-g_pituus_y * dy, g_pituus_y * dy, 2*g_pituus_y + 1);
    if numel(options.FWHM) == 3 && options.FWHM(3) > 0
        g_z = zeros(1,1,g_pituus_z*2+1);
        g_z(1,1,:) = linspace(-g_pituus_z * dz, g_pituus_z * dz, 2*g_pituus_z + 1);
        gaussK = gaussianKernel(g_x, g_y, g_z, options.FWHM(1) / (2 * sqrt(2 * log(2))), options.FWHM(2) / (2 * sqrt(2 * log(2))), options.FWHM(3) / (2 * sqrt(2 * log(2))));
    else
        gaussK = gaussianKernel(g_x, g_y, options.FWHM(1) / (2 * sqrt(2 * log(2))), options.FWHM(2) / (2 * sqrt(2 * log(2))));
    end
    
    if nargin >= 7
        options = uint32(g_pituus_x);
        varargout{1} = uint32(g_pituus_y);
        if numel(options.FWHM) == 3 && options.FWHM(3) > 0
            varargout{2} = uint32(g_pituus_z);
        end
    elseif nargin >= 3
        options = uint32(g_pituus_x);
        varargout{1} = uint32(g_pituus_y);
    else
        options.g_dim_x = uint32(g_pituus_x);
        options.g_dim_y = uint32(g_pituus_y);
        options.g_dim_z = uint32(g_pituus_z);
    end
    if isfield(options,'implementation') && (options.implementation == 2 || options.implementation == 3 || options.implementation == 5)
        gaussK = single(gaussK(:));
    end
else
    gaussK = single(0);
end
end

