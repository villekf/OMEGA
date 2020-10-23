function grad = AD(im, FluxType, Nx, Ny, Nz, options)
%AD Anisotropic Diffusion prior (AD) MRP
% Requires Image Processing Toolbox. Computes the MRP prior, but replaces
% the median filtered image with anisotropic diffusion smoothed image.
%
% Example:
%   grad = AD(im, FluxType, Nx, Ny, Nz, options)
% INPUTS:
%   im = The current estimate
%   FluxType = Flux type, 'quadratic' or 'exponential'
%   Nx = Image (estimate) size in X-direction
%   Ny = Image (estimate) size in Y-direction
%   Nz = Image (estimate) size in Z-direction
%   options = Gradient threshold for AD (KAD), Iteration count for AD (NiterAD)
%   options.med_no_norm = if true, then no normalization will be performed
%   (division by the AD smoothed image)
%
% OUTPUTS:
%   grad = The (gradient of) Anisotropic Diffusion smoothing prior
% 
% See also MRP, FMH, L_filter, Quadratic_prior, TVpriorfinal, TGV,
% Weighted_mean

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
if license('test', 'image_toolbox')
    grad = imdiffusefilt(reshape(im,Nx,Ny,Nz),'GradientThreshold', options.KAD, 'NumberOfIterations', options.NiterAD, 'ConductionMethod', FluxType);
    grad = grad(:);
    if options.med_no_norm
        grad = (im - grad);
    else
        grad = (im - grad) ./ grad;
    end
else
    warning('Image Processing Toolbox not found! Anisotropic diffusion can only be computed with Image Processing Toolbox license with implementations 1 or 4.')
end