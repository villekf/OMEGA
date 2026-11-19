function grad = Huber_prior(im, weights_huber, Nx, Ny, Nz, Ndx, Ndy, Ndz, delta)
%HUBER_PRIOR Huber Prior (HP)
%
% Example:
%   grad = Huber_prior(im, weights_huber, Nx, Ny, Nz, Ndx, Ndy, Ndz, delta)
% INPUTS:
%   im = The current estimate
%   weights_huber = Weights for huber prior
%   Nx = Image (estimate) size in X-direction
%   Ny = Image (estimate) size in Y-direction
%   Nz = Image (estimate) size in Z-direction
%   Ndx = How many pixels included in X-direction to weighting
%   Ndy = How many pixels included in Y-direction to weighting
%   Ndz = How many pixels included in Z-direction to weighting
%   delta = The Huber function delta parameter (constraint)
%
% OUTPUTS:
% grad = The Huber prior
%
% See also MRP, TGV, L_filter, FMH, TVpriorfinal, Weighted_mean,
% Quadratic_prior

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020  Ville-Veikko Wettenhovi, Samuli Summala
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
% The Huber gradient is grad_j = sum_k w_k * clip(x_j - x_k, -delta, delta), where
% k runs over the neighbors of voxel j. The weight kernel has the sum of the
% weights in the center and -w_k elsewhere, so w_k = -kernel(k). Each
% difference has to be clipped separately (clipping the weighted sum of the
% differences is not the Huber gradient). Matches the C++ Huber_prior.
im3 = reshape(im, Nx, Ny, Nz);
is2D = (Nz == 1 || Ndz == 0);
if Nz == 1
    ndzUsed = 0;
else
    ndzUsed = Ndz;
end
imPad = padding(im3, [Ndx Ndy ndzUsed]);
kx = 2 * Ndx + 1;
ky = 2 * Ndy + 1;
% A 2D image uses only the central axial slice of the kernel
if is2D
    kRange = Ndz;
else
    kRange = 0 : 2 * Ndz;
end
grad = zeros(Nx, Ny, Nz, 'like', im3);
maxAbs = 0;
for k = kRange
    for j = 0 : ky - 1
        for i = 0 : kx - 1
            if i == Ndx && j == Ndy && k == Ndz
                continue
            end
            w = -weights_huber(i + 1, j + 1, k + 1);
            if w == 0
                continue
            end
            if is2D
                shifted = imPad(i + 1 : i + Nx, j + 1 : j + Ny, :);
            else
                shifted = imPad(i + 1 : i + Nx, j + 1 : j + Ny, k + 1 : k + Nz);
            end
            diffIm = im3 - shifted;
            maxAbs = max(maxAbs, max(abs(diffIm(:))));
            grad = grad + w * min(max(diffIm, -delta), delta);
        end
    end
end
if maxAbs > 0 && delta >= maxAbs
    warning('Delta value of Huber prior larger than all the pixel difference values')
end
grad = grad(:);
