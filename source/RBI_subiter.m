function im = RBI_subiter(im, A, uu, epps, Summ, beta, dU, D, is_transposed)
%RBI_SUBITER Computes the Rescaled Block Iterative (RBI) estimate for the
%current subset/sub-iteration
%
% Example:
%   im = RBI_subiter(im, A, uu, epps, Summ, beta, dU, D)
% INPUTS:
%   im = The current estimate
%   A = The transpose of the (sparse) system matrix at current subset
%   uu = Measurements at current subset
%   epps = Small constant to prevent division by zero
%   Summ = sum(A,2)
%   beta = Regularization parameter (0, if no regularization)
%   dU = Gradient of the prior
%   D = 
%   is_transposed = true if A matrix is the transpose of it, false if not
%
% OUTPUTS:
%   im = The updated estimate
%
% See also BSREM_iter, BSREM_subiter, MBSREM, COSEM_OSL, ROSEM_subiter,
% OSL_OSEM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi
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
%%%% This version is from []
Summ = max((Summ + beta .* dU) ./ (D + beta .* dU));
if is_transposed
    im = im + (1/Summ).*(im./(D + beta .* dU)).*(A*(uu./(A'*im + epps) - 1) - beta .* dU);
else
    im = im + (1/Summ).*(im./(D + beta .* dU)).*(A'*(uu./(A*im + epps) - 1) - beta .* dU);
end
%%%% This is the original version
% Summ = max(Summ + beta .* dU);
% im = im + (im./Summ).*(A*(uu./(A'*im + epps) - 1) - beta .* dU);