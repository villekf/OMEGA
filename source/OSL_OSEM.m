function im = OSL_OSEM(im, Summ, beta, dU, epps, varargin)
%OSL_OSEM Computes the One-Step Late OSEM (OSL-OSEM) MAP estimates
%
% Example:
%   im = OSL_OSEM(im, Summ, beta, dU, A, uu, epps)
% INPUTS:
%   im = The current estimate
%   Summ = sum(A,2)
%   beta = Regularization parameter
%   dU = Gradient of the prior
%   epps = Small constant to prevent division by zero
%   A = The transpose of the (sparse) system matrix at current subset
%   uu = Measurements at current subset
%   is_transposed = true if A matrix is the transpose of it, false if not
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
if nargin == 8
    if varargin{3}
        im = (im./(Summ + beta * dU)).*(varargin{1} * (varargin{2}./(varargin{1}'*im + epps)));
    else
        im = (im./(Summ + beta * dU)).*(varargin{1}' * (varargin{2}./(varargin{1}*im + epps)));
    end
elseif nargin == 6
    im = (im./(Summ + beta * dU)).*(varargin{1} + epps);
else
    error('Invalid number of input arguments')
end