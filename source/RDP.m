function [grad, varargout] = RDP(im, weights_RDP, gamma, Nx, Ny, Nz, Ndx, Ndy, Ndz, tr_offsets)
%RDP Relative Difference Prior (RDP)
%   This function computes the gradient of the relative difference prior.
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
%   tr_offsets = Needed for vectorized indexing
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
im_orig = im(:);
if Nz==1
    im = padding(reshape(im,Nx,Ny),[Ndx Ndy]);
else
    im = padding(reshape(im,Nx,Ny,Nz),[Ndx Ndy Ndz]);
end
if nargout >= 2
    im = im(:);
    %     options = varargin{1};
    weights = [weights_RDP(1:ceil(numel(weights_RDP(:))/2) - 1),weights_RDP(ceil(numel(weights_RDP(:))/2) + 1 : end)]';
    apu = bsxfun(@minus,im(tr_offsets(:,ceil(numel(weights_RDP(:))/2))),im(tr_offsets(:,[1:(ceil(numel(weights_RDP(:))/2)-1) (ceil(numel(weights_RDP(:))/2)+1):end])))*weights;
    apu3 = bsxfun(@plus,im(tr_offsets(:,ceil(numel(weights_RDP(:))/2))),im(tr_offsets(:,[1:(ceil(numel(weights_RDP(:))/2)-1) (ceil(numel(weights_RDP(:))/2)+1):end])))*weights;
    mrf = -apu.^2 / (apu3 + gamma * abs(apu));
    varargout{1} = mrf;
end
try
    weights = [weights_RDP(1:ceil(numel(weights_RDP(:))/2) - 1);weights_RDP(ceil(numel(weights_RDP(:))/2) + 1 : end)];
    rjk = bsxfun(@rdivide,im(tr_offsets(:,ceil(numel(weights_RDP(:))/2))),im(tr_offsets(:,[1:(ceil(numel(weights_RDP(:))/2)-1) (ceil(numel(weights_RDP(:))/2)+1):end])));
%     apu = bsxfun(@minus,im(tr_offsets(:,ceil(numel(weights_RDP(:))/2))),im(tr_offsets(:,[1:(ceil(numel(weights_RDP(:))/2)-1) (ceil(numel(weights_RDP(:))/2)+1):end])));
%     apu2 = bsxfun(@plus,im(tr_offsets(:,ceil(numel(weights_RDP(:))/2))), 3 * im(tr_offsets(:,[1:(ceil(numel(weights_RDP(:))/2)-1) (ceil(numel(weights_RDP(:))/2)+1):end])));
%     apu3 = bsxfun(@plus,im(tr_offsets(:,ceil(numel(weights_RDP(:))/2))),im(tr_offsets(:,[1:(ceil(numel(weights_RDP(:))/2)-1) (ceil(numel(weights_RDP(:))/2)+1):end])));
%     grad = -((apu .* (gamma .* abs(apu) + apu2)) ./ (apu3 + gamma .* abs(apu)).^2) * weights;
%     grad = (((rjk - 1) .* (gamma .* abs(rjk - 1) + rjk + 3)) ./ (rjk + 1 + gamma .* abs(rjk - 1)).^2) * weights;
    grad = zeros(size(im_orig));
catch ME
    warning('Out of memory, switching to slower fallback code')
    weights_RDP(isinf(weights_RDP)) = 0;
    weights_RDP = reshape(weights_RDP, Ndy * 2 + 1, Ndx * 2 + 1, Ndz * 2 + 1);
    im_orig = reshape(im_orig,Nx,Ny,Nz);
    grad = zeros(size(im_orig));
    for kk = 1 : numel(grad)
        z = floor((kk - 1) / (Nx * Ny)) + 1;
        x = floor((kk - 1) / Nx) + 1 - (Nx * (z - 1));
        y = mod(kk - 1, Ny) + 1;
        for zz = -Ndz : Ndz
            for yy = -Ndy : Ndy
                for xx = -Ndx : Ndx
                    apu = im_orig(x, y, z) / im(x + xx + 1, y + yy + 1, z + zz + 1);
                    grad(x, y, z) = grad(x, y, z) + weights_RDP(xx + Ndx + 1,yy + Ndy + 1,zz + Ndz + 1) * ((apu - 1)*(gamma * abs(apu - 1) + apu + 3)) / (apu + 1 + gamma * abs(apu - 1)).^2;
                    %                     apu = im_orig(y, x, z) - im(y + yy + 1, x + xx + 1, z + zz + 1);
                    %                     grad(y, x, z) = -weights_RDP(yy + Ndy + 1,xx + Ndx + 1,zz + Ndz + 1) * (apu * (gamma * abs(apu) + im_orig(y, x, z) + 3 * im(y + yy + 1, x + xx + 1, z + zz + 1))) ...
                    %                         ./ (im_orig(y, x, z) + im(y + yy + 1, x + xx + 1, z + zz + 1) + gamma * abs(apu)).^2;
                end
            end
        end
    end
end
grad = grad(:);