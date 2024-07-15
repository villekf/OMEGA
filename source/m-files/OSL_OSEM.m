function im = OSL_OSEM(im, Summ, beta, dU, epps, varargin)
%OSL_OSEM Computes the One-Step Late OSEM (OSL-OSEM) MAP estimates
%
% Examples:
%   im = OSL_OSEM(im, Summ, beta, dU, epps, A, uu, SinDelayed, is_transposed)
%   im = OSL_OSEM(im, Summ, beta, dU, epps, RHS)
% INPUTS:
%   im = The current estimate
%   Summ = sum(A,2)
%   beta = Regularization parameter
%   dU = Gradient of the prior
%   epps = Small constant to prevent division by zero
%   A = The transpose of the (sparse) system matrix at current subset
%   uu = Measurements at current subset
%   SinDelayed = Randoms and/or scatter correction data. Dimension must be
%   either a scalar or a vector of same size as uu. If no scatter and/or
%   randoms data is available, use zero. 
%   is_transposed = true if A matrix is the transpose of it, false if not
%   RHS = The right hand side of OSEM (RHS = A'*(uu./(A*im + SinDelayed))) 
%
% OUTPUTS:
%   im = The updated estimate
%
% See also BSREM_iter, BSREM_subiter, MBSREM, COSEM_OSL, ROSEM_subiter,
% RBI_subiter

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
if nargin == 9
    if varargin{4}
        im = (im./(Summ + beta * dU + epps)).*(varargin{1} * (varargin{2}./(varargin{1}'*im + epps + varargin{3})) + epps);
    else
        im = (im./(Summ + beta * dU + epps)).*(varargin{1}' * (varargin{2}./(varargin{1}*im + epps + varargin{3})) + epps);
    end
elseif nargin == 6
%     temp = (Summ + beta * dU + epps);
%     temp(temp < 0) = epps;
    im = (im./(Summ + beta * dU + epps)).*(varargin{1});
%     joku = reshape(temp, 128, 128, 128);
%     joku2 = reshape(im./(temp), 128, 128, 128);
%     joku3 = reshape(im, 128, 128, 128);
%     joku4 = reshape(dU, 128, 128, 128);
%     subplot 221
%     imagesc(joku(:,:,64))
%     axis image
%     subplot 222
%     imagesc(joku2(:,:,64))
%     axis image
%     subplot 223
%     imagesc(joku3(:,:,64))
%     axis image
%     subplot 224
%     imagesc(joku4(:,:,64))
%     axis image
else
    error('Invalid number of input arguments')
end