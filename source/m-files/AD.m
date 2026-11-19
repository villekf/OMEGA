function grad = AD(im, FluxType, Nx, Ny, Nz, options)
%AD Anisotropic Diffusion prior (AD) MRP
% Computes the MRP prior, but replaces the median filtered image with
% anisotropic diffusion smoothed image. The smoothing is the same as
% af::anisotropicDiffusion of ArrayFire (used by the C++ implementation):
% 2D diffusion of each transaxial slice (first two dimensions), with the
% conductance normalized by the mean squared image gradient, replicated
% borders and simultaneous (Jacobi) updates of all voxels as in the
% OpenCL/CUDA backends of ArrayFire.
%
% Example:
%   grad = AD(im, FluxType, Nx, Ny, Nz, options)
% INPUTS:
%   im = The current estimate
%   FluxType = Flux type, 'quadratic' (or 2) or 'exponential' (or 1)
%   Nx = Image (estimate) size in X-direction
%   Ny = Image (estimate) size in Y-direction
%   Nz = Image (estimate) size in Z-direction
%   options = Gradient threshold for AD (KAD), Iteration count for AD
%   (NiterAD), time step (TimeStepAD), diffusion type (DiffusionType, 1 =
%   gradient-based, 2 = modified curvature diffusion equation (MCDE),
%   which always uses the exponential flux), epps
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
if ischar(FluxType) || isstring(FluxType)
    quadFlux = strcmpi(FluxType, 'quadratic');
else
    quadFlux = FluxType == 2;
end
MCDE = isfield(options, 'DiffusionType') && options.DiffusionType == 2;
x = reshape(im, Nx, Ny, Nz);
dt = double(options.TimeStepAD);
K = double(options.KAD);
cnst = -2 * K * K / numel(x);
ip = [2 : Nx, Nx];
im1 = [1, 1 : Nx - 1];
jp = [2 : Ny, Ny];
jm = [1, 1 : Ny - 1];
for it = 1 : double(options.NiterAD)
    % The conductance is normalized with the mean squared gradient magnitude (central differences, one-sided at the borders)
    g0 = gradAF(x, 1);
    g1 = gradAF(x, 2);
    mct = 1 / (cnst * sum(g0(:) .* g0(:) + g1(:) .* g1(:)));
    C = x;
    E = x(ip, :, :);
    W = x(im1, :, :);
    S = x(:, jp, :);
    N = x(:, jm, :);
    SE = x(ip, jp, :);
    NE = x(ip, jm, :);
    SW = x(im1, jp, :);
    NW = x(im1, jm, :);
    dx = (E - W) * 0.5;
    dy = (S - N) * 0.5;
    % First dimension
    df0 = E - C;
    db0 = C - W;
    gmsqf = df0 .* df0 + 0.25 * (dy + 0.5 * (SE - NE)).^2;
    gmsqb = db0 .* db0 + 0.25 * (dy + 0.5 * (SW - NW)).^2;
    % Second dimension
    df = S - C;
    db = C - N;
    gmsqf2 = df .* df + 0.25 * (dx + 0.5 * (SE - SW)).^2;
    gmsqb2 = db .* db + 0.25 * (dx + 0.5 * (NE - NW)).^2;
    if MCDE
        delta = (df0 ./ sqrt(1e-10 + gmsqf)) .* exp(gmsqf * mct) - (db0 ./ sqrt(1e-10 + gmsqb)) .* exp(gmsqb * mct);
        delta = delta + (df ./ sqrt(1e-10 + gmsqf2)) .* exp(gmsqf2 * mct) - (db ./ sqrt(1e-10 + gmsqb2)) .* exp(gmsqb2 * mct);
        pos = delta > 0;
        propGrad = min(db0, 0).^2 + max(df0, 0).^2 + min(db, 0).^2 + max(df, 0).^2;
        propGradN = max(db0, 0).^2 + min(df0, 0).^2 + max(db, 0).^2 + min(df, 0).^2;
        propGrad(~pos) = propGradN(~pos);
        delta = sqrt(propGrad) .* delta;
    else
        if quadFlux
            delta = (1 ./ (1 + gmsqf * mct)) .* df0 - (1 ./ (1 + gmsqb * mct)) .* db0;
            delta = delta + (1 ./ (1 + gmsqf2 * mct)) .* df - (1 ./ (1 + gmsqb2 * mct)) .* db;
        else
            delta = exp(gmsqf * mct) .* df0 - exp(gmsqb * mct) .* db0;
            delta = delta + exp(gmsqf2 * mct) .* df - exp(gmsqb2 * mct) .* db;
        end
    end
    x = C + delta * dt;
end
grad = x(:);
if options.med_no_norm
    grad = (im - grad);
else
    grad = (im - grad) ./ (grad + options.epps);
end
end

function g = gradAF(x, dim)
% Gradient along the first (dim = 1) or second (dim = 2) dimension as in af::grad: central differences inside, one-sided
% differences at the borders
n = size(x, dim);
g = zeros(size(x), 'like', x);
if n < 2
    return
end
if dim == 1
    g(2 : n - 1, :, :) = 0.5 * (x(3 : n, :, :) - x(1 : n - 2, :, :));
    g(1, :, :) = x(2, :, :) - x(1, :, :);
    g(n, :, :) = x(n, :, :) - x(n - 1, :, :);
else
    g(:, 2 : n - 1, :) = 0.5 * (x(:, 3 : n, :) - x(:, 1 : n - 2, :));
    g(:, 1, :) = x(:, 2, :) - x(:, 1, :);
    g(:, n, :) = x(:, n, :) - x(:, n - 1, :);
end
end
