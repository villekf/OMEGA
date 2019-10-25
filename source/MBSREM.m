function im = MBSREM(im, varargin)
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
% Copyright (C) 2019  Ville-Veikko Wettenhovi, Samuli Summala
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
U = varargin{1};
epps = varargin{4};
epsilon_mramla = varargin{6};
iter = varargin{8};
randoms_correction = varargin{10};
is_transposed = varargin{11};

pp = im < (U/2);
UU = zeros(size(im,1),1);
UU(pp) = im(pp)./(varargin{2}(pp)+epps);
UU(~pp) = (U-im(~pp))./(varargin{2}(~pp)+epps);
if is_transposed
    if randoms_correction
        l_mramla = varargin{3}'*im + varargin{9};
    else
        l_mramla = varargin{3}'*im;
    end
else
    if randoms_correction
        l_mramla = varargin{3}*im + varargin{9};
    else
        l_mramla = varargin{3}*im;
    end
end
l_mramla(l_mramla <= 0) = epps;
if randoms_correction
    pp = l_mramla <= epsilon_mramla & varargin{5} > 0 & varargin{9} == 0;
else
    pp = l_mramla <= epsilon_mramla & varargin{5} > 0;
end
hr = varargin{5}./l_mramla - 1;
hr(pp) = varargin{5}(pp)./epsilon_mramla - 1 - (varargin{5}(pp)./epsilon_mramla.^2).*(l_mramla(pp) - epsilon_mramla);
if is_transposed
    if nargin > 12
        im = im + varargin{7}(iter).*UU.*(varargin{3}*hr - varargin{12} .* varargin{13} + epps);
    else
        im = im + varargin{7}(iter).*UU.*(varargin{3}*hr + epps);
    end
else
    if nargin > 12
        im = im + varargin{7}(iter).*UU.*(varargin{3}'*hr - varargin{12} .* varargin{13} + epps);
    else
        im = im + varargin{7}(iter).*UU.*(varargin{3}'*hr + epps);
    end
end
im(im <= 0) = epps;
im(im >= U) = U - epps;