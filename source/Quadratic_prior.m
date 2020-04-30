function grad = Quadratic_prior(im, ~, ~, weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz)
%QUADRATIC_PRIOR Quadratic Prior (QP)
%
% Example:
%   grad = Quadratic_prior(im, tr_offsets, weights, weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz)
% INPUTS:
%   im = The current estimate
%   tr_offsets = Offset values, based on Ndx, Ndy and Ndz
%   weights = Distance-based weights
%   weights_quad = Weights for quadratic prior
%   Nx = Image (estimate) size in X-direction
%   Ny = Image (estimate) size in Y-direction
%   Nz = Image (estimate) size in Z-direction
%   Ndx = How many pixels included in X-direction to weighting
%   Ndy = How many pixels included in Y-direction to weighting
%   Ndz = How many pixels included in Z-direction to weighting
%
% OUTPUTS:
% grad = The quadratic prior
%
% See also MRP, TGV, L_filter, FMH, TVpriorfinal, Weighted_mean

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ville-Veikko Wettenhovi, Samuli Summala
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
if Nz==1
    im = padding(reshape(im,Nx,Ny),[Ndx Ndy]);
else
    im = padding(reshape(im,Nx,Ny,Nz),[Ndx Ndy Ndz]);
end
grad = convn(im, weights_quad, 'valid');
grad = grad(:);
% grad = bsxfun(@minus,im(tr_offsets(:,isinf(weights))),im(tr_offsets(:,[1:(find(isinf(weights))-1) (find(isinf(weights))+1):end])))*weights_quad;