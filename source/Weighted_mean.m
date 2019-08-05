function grad = Weighted_mean(im, tr_offsets, weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, mean_type, epps, w_sum, med_no_norm)
%Weighted_mean Weighted mean prior
%
% Example:
%   grad = Weighted_mean(im, tr_offsets, weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, mean_type, epps, w_sum, med_no_norm)
% INPUTS:
%   im = The current estimate
%   tr_offsets = Offset values, based on Ndx, Ndy and Ndz
%   weighted_weights = Weights for the weighted mean
%   Nx = Image (estimate) size in X-direction
%   Ny = Image (estimate) size in Y-direction
%   Nz = Image (estimate) size in Z-direction
%   Ndx = How many pixels included in X-direction to weighting
%   Ndy = How many pixels included in Y-direction to weighting
%   Ndz = How many pixels included in Z-direction to weighting
%   mean_type = Type of the mean, 1 = Arithmetic, 2 = Harmonic, 3 = Geometric
%   epps = Small constant to prevent division by zero
%   w_sum = Sum of the weighted_weights
%   med_no_norm = if true, then no normalization will be performed
%   (division by the weighted mean filtered image)
%
% OUTPUTS:
%   grad = The (gradient of) weighted mean prior
%
% See also MRP, TGV, L_filter, Quadratic_prior, TVpriorfinal, FMH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi, Samuli Summala
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
    padd = padding(reshape(im,Nx,Ny),[Ndx Ndy]);
else
    padd = padding(reshape(im,Nx,Ny,Nz),[Ndx Ndy Ndz]);
end
if mean_type == 1
    grad = (padd(tr_offsets)*weighted_weights)./w_sum;
elseif mean_type == 2
    grad = w_sum./((1./padd(tr_offsets))*weighted_weights);
elseif mean_type == 3
    grad = exp((log(padd(tr_offsets))*weighted_weights)./w_sum);
else
    error('Unsupported mean')
end
if med_no_norm
    grad = (im - grad);
else
    grad = (im - grad)./(grad + epps);
end