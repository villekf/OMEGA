function output = MBSREM(im, rhs, options, U, lam, iter, osa_iter, ii)
%MBSREM Computes the modified Block Sequential Regularized Expectation
% Maximization (MBSREM) estimates
%
% Examples:
%   im = MBSREM(im, U, D, A, epps, uu, epsilon_mramla, lambda, iter, SinDelayed, randoms_correction, is_transposed)
%   im = MBSREM(im, U, D, A, epps, uu, epsilon_mramla, lambda, iter, SinDelayed, randoms_correction, is_transposed, beta, dU)
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
%   SinDelayed = Randoms and/or scatter correction data. Dimension must be
%   either a scalar or a vector of same size as varargin{5}. If no scatter
%   and/or randoms data is available, use zero.
%   randoms_correction = Boolean parameter that is true when randoms and/or
%   scatter correction is applied and false otherwise.
%   is_transposed = true if varargin{3} matrix is the transpose of it,
%   false if not
%   beta = Regularization parameter (0, if no regularization)
%   dU = Gradient of the prior
%
% OUTPUTS:
%   im = The updated estimate
%
% See also BSREM_iter, BSREM_subiter, MBSREM_epsilon, MBSREM_prepass,
% RBI_subiter, COSEM_OSL, ROSEM_subiter, OSL_OSEM

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
% MERCHANTABILITY or FITNESS FOR varargin{3} PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kk = (iter - 1) * options.subsets + osa_iter;
pp = im >= (U / 2);
if (any(pp))
    apuIm = im;
    apuIm(pp) = U - apuIm(pp);
    rhs = applyImagePreconditioning(options, rhs, apuIm, kk, ii);
else
    rhs = applyImagePreconditioning(options, rhs, im, kk, ii);
end
output = im + lam(iter) .* rhs;
output(output < options.epps) = options.epps;
output(output >= U) = U - options.epps;
end