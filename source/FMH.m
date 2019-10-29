function grad = FMH(im, tr_offsets, fmh_weights, weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, med_no_norm)
%FMH Finite Impulse Response (FIR) Median Hybrid prior (FMH)
%
% Example:
%   grad = FMH(im, tr_offsets, fmh_weights, weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, med_no_norm)
% INPUTS:
%   im = The current estimate
%   tr_offsets = Offset values, based on Ndx, Ndy and Ndz
%   fmh_weights = Weights for FMH
%   Nx = Image (estimate) size in X-direction
%   Ny = Image (estimate) size in Y-direction
%   Nz = Image (estimate) size in Z-direction
%   Ndx = How many pixels included in X-direction to weighting
%   Ndy = How many pixels included in Y-direction to weighting
%   Ndz = How many pixels included in Z-direction to weighting
%   epps = Small constant to prevent division by zero
%   med_no_norm = if true, then no normalization will be performed
%   (division by the FMH filtered image)
%
% OUTPUTS:
%   grad = The (gradient of) FMH prior
%
% See also MRP, TGV, L_filter, Quadratic_prior, TVpriorfinal, Weighted_mean

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
alku = (find(isinf(weights)) - 1);
if Nz == 1 || Ndz == 0
    grad = zeros(N, 5);
    luup = 4;
else
    grad = zeros(N, 14);
    luup = 13;
end
for ii = 1 : luup
    grad(:,ii) = padd(tr_offsets(:,1+Ndx*(ii-1):alku/Ndx-(ii-1):end - Ndx*(ii-1)))*fmh_weights(:,ii);
end
grad(:,end) = padd(tr_offsets(:,isinf(weights)));
grad = median(grad,2);
if med_no_norm
    grad = (im - grad);
else
    grad = (im - grad)./(grad + epps);
end