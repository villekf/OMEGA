function grad = L_filter(im, tr_offsets, a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, varargin)
%L_FILTER L-filter prior (L-filter)
%   Computes the "gradient" of the L-filter. Behaves identically to MRP,
%   except that median has been replaced by the L-filter.
%
% Examples:
%   grad = L_filter(im, tr_offsets, a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz)
%   grad = L_filter(im, tr_offsets, a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, med_no_norm)
% INPUTS:
%   im = The current estimate
%   tr_offsets = Offset values, based on Ndx, Ndy and Ndz, computed by
%   computeOffsets
%   a_L = Weights for L-filter (Laplace distributed)
%   Nx = Image (estimate) size in X-direction
%   Ny = Image (estimate) size in Y-direction
%   Nz = Image (estimate) size in Z-direction
%   Ndx = How many pixels included in X-direction to weighting
%   Ndy = How many pixels included in Y-direction to weighting
%   Ndz = How many pixels included in Z-direction to weighting
%   epps = Small constant to prevent division by zero (optional). Default
%   value is 1e-8.
%   med_no_norm = if true, then no normalization will be performed
%   (division by the L-filtered image) (optional). Default value is false.
%
% OUTPUTS:
% grad = The (gradient of) L-filter
%
% See also MRP, TGV, FMH, Quadratic_prior, TVpriorfinal, Weighted_mean,
% lfilter_weights, computeOffsets

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
if nargin >= 10 && ~isempty(varargin{1})
    epps = varargin{1};
else
    epps = 1e-8;
end
if nargin >= 11 && ~isempty(varargin{2})
    med_no_norm = varargin{2};
else
    med_no_norm = false;
end
if Nz==1
    padd = padding(reshape(im,Nx,Ny),[Ndx Ndy]);
else
    padd = padding(reshape(im,Nx,Ny,Nz),[Ndx Ndy Ndz]);
end
grad = sort(padd(tr_offsets),2)*a_L;
if med_no_norm % Do not compute normalization
    grad = (im - grad);
else % Compute normalization
    grad = (im - grad)./(grad + epps);
end