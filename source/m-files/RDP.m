function grad = RDP(im, ~, gamma, Nx, Ny, Nz, ~, ~, ~, ~, varargin)
%RDP Relative Difference Prior (RDP)
%   This function computes the gradient of the relative difference prior.
%   The standard RDP uses only the 6 face neighbors with unit weights and
%   clamp-to-edge (Neumann) boundary conditions, exactly as the C++/GPU
%   implementation (non-corner RDP). The weights, Ndx, Ndy, Ndz and
%   tr_offsets inputs are not used, they are kept for compatibility.
%
% Example:
%   grad = RDP(im, weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz)
% INPUTS:
%   im = The current estimate
%   weights_RDP = Weights for RDP
%   gamma = The edge preservation coefficient
%   Nx = Image (estimate) size in X-direction
%   Ny = Image (estimate) size in Y-direction
%   Nz = Image (estimate) size in Z-direction
%   Ndx = How many pixels included in X-direction neighborhood
%   Ndy = How many pixels included in Y-direction neighborhood
%   Ndz = How many pixels included in Z-direction neighborhood
%   tr_offsets = Not used
%   epps = Small value used in the denominator (optional, default 1e-6)
%
% OUTPUTS:
% grad = The relative difference prior
%
% See also MRP, TGV, L_filter, FMH, TVpriorfinal, Weighted_mean,
% Quadratic_prior

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
if nargin >= 11 && ~isempty(varargin{1})
    epps = varargin{1};
else
    epps = 1e-6;
end
im = reshape(im, Nx, Ny, Nz);
epps = cast(epps, 'like', im);
gamma = cast(gamma, 'like', im);
% Neighbor indices clamped to the volume (clamp-to-edge)
xp = min((1:Nx) + 1, Nx);
xm = max((1:Nx) - 1, 1);
yp = min((1:Ny) + 1, Ny);
ym = max((1:Ny) - 1, 1);
zp = min((1:Nz) + 1, Nz);
zm = max((1:Nz) - 1, 1);
% Positive (x + 1, y + 1, z + 1) and negative (x - 1, y - 1, z - 1) side neighbors. NaN is checked separately
% for the sums of the positive and negative sides, as in the C++ code
outPos = rdpTerm(im, im(xp, :, :), gamma, epps) + rdpTerm(im, im(:, yp, :), gamma, epps) + rdpTerm(im, im(:, :, zp), gamma, epps);
outNeg = rdpTerm(im, im(xm, :, :), gamma, epps) + rdpTerm(im, im(:, ym, :), gamma, epps) + rdpTerm(im, im(:, :, zm), gamma, epps);
outPos(isnan(outPos)) = 0;
outNeg(isnan(outNeg)) = 0;
grad = outPos(:) + outNeg(:);
end

function out = rdpTerm(uj, uk, gamma, epps)
% One neighbor term of the RDP gradient
d = uj - uk;
divPow2 = uj + uk + gamma .* abs(d);
out = d .* (gamma .* abs(d) + uj + 3 .* uk + epps .* epps) ./ (divPow2 .* divPow2 + epps);
end
