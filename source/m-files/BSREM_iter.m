function im = BSREM_iter(im, lambda, iter, beta, dU, epps)
%BSREM_ITER Computes the Block Sequential Regularized Expectation Maximization
% (BSREM) or ROSEM-MAP estimates for current iteration
%
% Example:
%   im = BSREM_iter(im, lambda, iter, beta, dU, epps)
% INPUTS:
%   im = The current estimate
%   lambda = Relaxation parameter
%   iter = Current iteration number
%   beta = Regularization parameter (0, if no regularization)
%   dU = Gradient of the prior
%   epps = Small constant to prevent negative values
%
% OUTPUTS:
%   im = The updated estimate
%
% See also BSREM_subiter, MBSREM, COSEM_OSL, ROSEM_subiter, RBI_subiter, 
% OSL_OSEM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019 Ville-Veikko Wettenhovi, Samuli Summala
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
im = im - beta.*lambda(iter).*im.*dU;
im(im < epps) = epps;