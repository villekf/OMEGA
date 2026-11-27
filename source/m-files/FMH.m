function grad = FMH(im, tr_offsets, fmh_weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, varargin)
%FMH Finite Impulse Response (FIR) Median Hybrid prior (FMH)
%   Computes the "gradient" of the FMH. Behaves identically to MRP,
%   except that median has been replaced by the FMH.
%
% Examples:
%   grad = FMH(im, tr_offsets, fmh_weights, weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz)
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
%   epps = Small constant to prevent division by zero (optional). Default
%   value is 1e-8.
%   med_no_norm = if true, then no normalization will be performed
%   (division by the FMH filtered image) (optional). Default value is
%   false.
%
% OUTPUTS:
%   grad = The (gradient of) FMH prior
%
% See also MRP, TGV, L_filter, Quadratic_prior, TVpriorfinal, Weighted_mean

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
if nargin >= 12 && ~isempty(varargin{1})
    epps = varargin{1};
else
    epps = 1e-8;
end
if nargin >= 13 && ~isempty(varargin{2})
    med_no_norm = varargin{2};
else
    med_no_norm = false;
end
if Nz==1
    padd = padding(reshape(im,Nx,Ny),[Ndx Ndy]);
else
    padd = padding(reshape(im,Nx,Ny,Nz),[Ndx Ndy Ndz]);
end
alku = floor(size(tr_offsets,2) / 2);
if Nz == 1 || Ndz == 0
    grad = zeros(N, 5);
    luup = 4;
else
    grad = zeros(N, 14);
    luup = 13;
end
is2DFMH = (Nz == 1 || Ndz == 0);
for ii = 1 : luup
    % Line direction ii in the unit cube (2D: 4 directions, 3D: 13 directions). The neighborhood is the
    % (2Ndx+1)x(2Ndy+1)x(2Ndz+1) box with x fastest, so the columns of the line are the center plus k * d, k = N...-N.
    % Requires Ndx == Ndy and Ndz equal to 0 or Ndx. Matches the C++ FMH.
    id0 = ii - 1;
    dx = mod(id0, 3) - 1;
    if is2DFMH
        dy = floor(id0 / 3) - 1;
        dz = 0;
    else
        dy = mod(floor(id0 / 3), 3) - 1;
        dz = floor(id0 / 9) - 1;
    end
    kk = double(Ndx) : -1 : -double(Ndx);
    colIdx = 1 + (double(Ndx) + kk * dx) + (double(Ndy) + kk * dy) * (2 * double(Ndx) + 1) ...
        + (double(Ndz) + kk * dz) * ((2 * double(Ndx) + 1) * (2 * double(Ndy) + 1));
    grad(:,ii) = padd(tr_offsets(:,colIdx))*fmh_weights(:,ii);
end
grad(:,end) = padd(tr_offsets(:,alku + 1));
grad = median(grad,2);
if med_no_norm
    grad = (im - grad);
else
    grad = (im - grad)./(grad + epps);
end