function sens = OSL(Summ, beta, dU, epps)
%OSL Computes the One-Step Late regularization
%
% Examples:
%   grad = OSL(Summ, beta, dU, epps, is_transposed)
% INPUTS:
%   Summ = sum(A,2)
%   beta = Regularization parameter
%   dU = Gradient of the prior
%   epps = Small constant to prevent negative values or zeroes
%
% OUTPUTS:
%   sens = The regularized sensitivity image
%
% See also BSREM_iter, BSREM_subiter, MBSREM, COSEM_OSL, ROSEM_subiter,
% RBI_subiter

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
if nargin == 4
    sens = (Summ + beta * dU + epps);
%     sens(sens < epps) = epps;
else
    error('Invalid number of input arguments')
end