function grad = MRP(im, medx, medy, medz, Nx, Ny, Nz, epps, tr_offsets, med_no_norm)
%MRP Median Root Prior (MRP)
%
% Example:
%   grad = MRP(im, medx, medy, medz, Nx, Ny, Nz, epps, tr_offsets, med_no_norm)
% INPUTS:
%   im = The current estimate
%   medx = Size of median in X-direction
%   medy = Size of median in Y-direction
%   medz = Size of median in Z-direction
%   Nx = Image (estimate) size in X-direction
%   Ny = Image (estimate) size in Y-direction
%   Nz = Image (estimate) size in Z-direction
%   epps = Small constant to prevent division by zero
%   tr_offsets = Offset values, based on Ndx, Ndy and Ndz (only used if
%   image processing toolbox is not present)
%   med_no_norm = if true, then no normalization will be performed
%   (division by the median filtered image)
%
% OUTPUTS:
% grad = The (gradient of) MRP
%
% See also FMH, TGV, L_filter, Quadratic_prior, TVpriorfinal, Weighted_mean

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
% If Image Processing Toolbox exists
if license('test', 'image_toolbox')
    if Nz==1
        grad = medfilt2(reshape(im,Nx,Ny), [medx medy], 'symmetric');
    else
        grad = medfilt3(reshape(im,Nx,Ny,Nz), [medx medx medz]);
    end
    % If Image Processing Toolbox does not exist or Octave is used (slower)
else
    if Nz==1
        padd = padding(reshape(im, Nx, Ny), [max(0, floor(medx / 2)) max(0, floor(medy / 2))]);
    else
        padd = padding(reshape(im, Nx, Ny, Nz), [max(0, floor(medx / 2)) max(0, floor(medy / 2)) max(0, floor(medz / 2))]);
    end
    grad = padd(tr_offsets);
    grad = median(grad, 2);
end
grad = grad(:);
if med_no_norm
    grad = (im - grad);
else
    grad = (im - grad) ./ (grad + epps);
end