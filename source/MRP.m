function grad = MRP(im, medx, medy, medz, Nx, Ny, Nz, varargin)
%MRP Median Root Prior (MRP)
%   Computes the "gradient" of the median root prior in the [medx medy
%   medz] neighorbood. Uses the MATLAB function medfilt3 if it is
%   available. If not, then uses the offsets values computed by
%   computeOffsets (tr_offsets) to compute the median. Can be used with or
%   without the normalization in MRP (see below or wiki for more
%   information).
%
% Examples:
%   grad = MRP(im, medx, medy, medz, Nx, Ny, Nz)
%   grad = MRP(im, medx, medy, medz, Nx, Ny, Nz, epps, tr_offsets, med_no_norm)
% INPUTS:
%   im = The current estimate
%   medx = Size of median region in X-direction (default medx = Ndx * 2 +
%   1) 
%   medy = Size of median region in Y-direction
%   medz = Size of median region in Z-direction
%   Nx = Image (estimate) size in X-direction
%   Ny = Image (estimate) size in Y-direction
%   Nz = Image (estimate) size in Z-direction
%   epps = Small constant to prevent division by zero (optional, default
%   value is 1e-8)
%   tr_offsets = Offset values, based on Ndx, Ndy and Ndz (only used if
%   image processing toolbox is not present, can be omitted if toolbox is
%   present)
%   med_no_norm = if true, then no normalization will be performed
%   (division by the median filtered image). Can be omitted, default is
%   false.
%
% OUTPUTS:
% grad = The (gradient of) MRP
%
% See also FMH, TGV, L_filter, Quadratic_prior, TVpriorfinal,
% Weighted_mean, computeOffsets

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
% If Image Processing Toolbox exists
if nargin >= 8 && ~isempty(varargin{1})
    epps = varargin{1};
else
    epps = 1e-8;
end
if nargin >= 10 && ~isempty(varargin{3})
    med_no_norm = varargin{3};
else
    med_no_norm = false;
end

if license('test', 'image_toolbox') || exist('medfilt3','file') == 2
    if Nz==1
        grad = medfilt2(reshape(im,Nx,Ny), [medx medy], 'symmetric');
    else
        grad = medfilt3(reshape(im,Nx,Ny,Nz), [medx medx medz]);
    end
else % If Image Processing Toolbox does not exist or Octave is used (slower and more memory intensive)
    if nargin >= 9 && isempty(varargin{2})
        error('tr_offsets variable needs to be included')
    end
    if Nz==1
        padd = padding(reshape(im, Nx, Ny), [max(0, floor(medx / 2)) max(0, floor(medy / 2))]);
    else
        padd = padding(reshape(im, Nx, Ny, Nz), [max(0, floor(medx / 2)) max(0, floor(medy / 2)) max(0, floor(medz / 2))]);
    end
    if min(varargin{2}(:) == 0)
        grad = padd(varargin{2} + uint32(1));
    else
        grad = padd(varargin{2});
    end
    grad = median(grad, 2);
end
grad = grad(:);
if med_no_norm % Do not compute normalization
    grad = (im - grad);
else % Compute normalization
    grad = (im - grad) ./ (grad + epps);
end