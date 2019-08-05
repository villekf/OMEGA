function im = MBSREM(im, U, D, A, epps, uu, epsilon_mramla, lambda, iter, beta, dU, is_transposed)
%MBSREM Computes the modified Block Sequential Regularized Expectation
% Maximization (MBSREM) estimates
%
% Example:
%   im = MBSREM(im, U, D, A, epps, uu, epsilon_mramla, lambda, iter, beta, dU)
% INPUTS:
%   im = The current estimate
%   U = Upper bound of MBSREM (see MBSREM_epsilon)
%   D = Sum of the complete data system matrix divided by the number of
%   subsets (D = sum(B,2)/subsets, where B = [A_1,A_2,...,A_subsets]
%   A = The transpose of the (sparse) system matrix at current subset
%   epps = Small constant to prevent division by zero
%   uu = Measurements at current subset
%   epsilon_mramla = Epsilon value of MBSREM (see MBSREM_epsilon)
%   lambda = Relaxation parameter
%   iter = Current iteration number
%   beta = Regularization parameter (0, if no regularization)
%   dU = Gradient of the prior
%
% OUTPUTS:
%   im = The updated estimate
%
% See also BSREM_iter, BSREM_subiter, MBSREM_epsilon, MBSREM_prepass,
% RBI_subiter, COSEM_OSL, ROSEM_subiter, OSL_OSEM

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
pp = im < (U/2);
UU = zeros(size(im,1),1);
UU(pp) = im(pp)./(D(pp)+epps);
UU(~pp) = (U-im(~pp))./(D(~pp)+epps);
if is_transposed
    l_mramla = A'*im;
else
    l_mramla = A*im;
end
l_mramla(l_mramla <= 0) = epps;
pp = l_mramla < epsilon_mramla;
hr = uu./l_mramla - 1;
hr(pp) = uu(pp)./l_mramla(pp) - 1 - (uu(pp)./l_mramla(pp).^2).*(l_mramla(pp) - epsilon_mramla);
if is_transposed
    im = im + lambda(iter).*UU.*(A*hr - beta .* dU);
else
    im = im + lambda(iter).*UU.*(A'*hr - beta .* dU);
end
im(im < 0) = epps;