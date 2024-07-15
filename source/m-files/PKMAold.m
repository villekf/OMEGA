function im = PKMA(im, options, A, Summ, lambda, alpha, sigma, epps, grad, beta, iter, osa_iter, varargin)
%PKMA Preconditioned Krasnoselskii-Mann algorithm with EM preconditioner
%   Computes the PKM algorithm with the EM preconditioner.
%
% Examples:
%   im = PKMA(im, options, A, Summ, lambda, alpha, sigma, epps, grad, beta, iter, osa_iter, uu, is_transposed, SinD, gaussK, Nx, Ny, Nz)
%   im = PKMA(im, options, RHS, Summ, lambda, alpha, sigma, epps, grad, beta, iter, sub_iter)
% INPUTS:
%   im = The current estimate
%   A = The transpose of the (sparse) system matrix at current subset
%   Summ = sum(A,2), sensitivity image
%   lambda = Relaxation parameter
%   alpha = Step size parameter
%   sigma = Additiona step size parameter
%   epps = Small constant to prevent division by zero
%   grad = Gradient of the prior
%   beta = Regularization parameter
%   iter = Current iteration number
%   sub_iter = Current subiteration
%   uu = Measurements at current subset
%   is_transposed = true if A matrix is the transpose of it, false if not
%   SinDelayed = Randoms and/or scatter correction data. Dimension must be
%   either a scalar or a vector of same size as varargin{5}. If no scatter
%   and/or randoms data is available, use zero. 
%   RHS = The backprojection
%
% OUTPUTS:
%   im = The updated estimate
%
% See also OSEM_im, DRAMA_subiter, ROSEM_subiter, RBI_subiter, COSEM_im,
% ACOSEM_im, ECOSEM_im, MBSREM, BSREM_subiter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021 Ville-Veikko Wettenhovi
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

S = (im + epps) ./ options.pj3;
im_ = im;
% Implementation 1
if nargin >= 13
    if options.use_psf
        im_apu = computeConvolution(im, options, varargin{4}, varargin{5}, varargin{6}, varargin{7});
    else
        im_apu = im;
    end
    if varargin{2}
        FP = (A' * im_apu + varargin{3});
        if options.CT
            FP = exp(FP);
        end
        BP = A * (ones(length(varargin{1}),1,'double') - varargin{1} ./ FP);
    else
        FP = (A * im + varargin{3});
        if options.CT
            FP = exp(FP);
        end
        BP = A' * (ones(length(varargin{1}),1,'double') - varargin{1} ./ FP);
    end
    if options.use_psf
        BP = computeConvolution(BP, options, varargin{4}, varargin{5}, varargin{6}, varargin{7});
    end
% Implementation 4
else
    BP = Summ(:,osa_iter) - A;
end
im = im - lambda(iter) .* S .* (BP + beta * grad);
im(im < epps) = epps;
ind = (iter - 1) * options.subsets + osa_iter;
im = (1 - alpha(ind)) .* im_ + alpha(ind) .* (sigma(ind) .* im);