function [grad, varargout] = Quadratic_prior(im, weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz, varargin)
%QUADRATIC_PRIOR Quadratic Prior (QP)
%   This function computes the gradient of the quadratic prior. The
%   MRF-function can be output as the (optional) second output.
%   MRF-function output requires tr_offsets.
%
% Example:
%   grad = Quadratic_prior(im, weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz)
% INPUTS:
%   im = The current estimate
%   weights_quad = Weights for quadratic prior
%   Nx = Image (estimate) size in X-direction
%   Ny = Image (estimate) size in Y-direction
%   Nz = Image (estimate) size in Z-direction
%   Ndx = How many pixels included in X-direction neighborhood
%   Ndy = How many pixels included in Y-direction neighborhood
%   Ndz = How many pixels included in Z-direction neighborhood
%   options = Only used when the number of output arguments are two. Needs
%   to contain tr_offsets-variable.
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
weights_quad = -weights_quad;
weights_quad(ceil(numel(weights_quad)/2)) = -weights_quad(ceil(numel(weights_quad)/2));
grad = convn(im, weights_quad, 'valid');
grad = grad(:);
if nargout >= 2
    im = im(:);
    options = varargin{1};
    weights_quad2 = abs(weights_quad(:));
    weights_quad2(ceil(numel(weights_quad2)/2)) = [];
    mrf = 0.5 .* bsxfun(@minus,im(options.tr_offsets(:,ceil(size(options.tr_offsets,2)/2))),im(options.tr_offsets(:,[1:(ceil(size(options.tr_offsets,2)/2)-1) (ceil(size(options.tr_offsets,2)/2)+1):end]))).^2*weights_quad2;
    varargout{1} = mrf;
end
% grad1 = bsxfun(@minus,im(tr_offsets(:,isinf(weights))),im(tr_offsets(:,[1:(find(isinf(weights))-1) (find(isinf(weights))+1):end])))*weights_quad;