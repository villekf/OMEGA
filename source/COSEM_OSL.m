function [im, C_co] = COSEM_OSL(im, D, beta, dU, A, uu, epps, C_co, h, COSEM_MAP, osa_iter, SinDelayed, is_transposed, varargin)
%COSEM_OSL Computes the One-Step Late COSEM (OSL-COSEM) MAP estimates
%
% Example:
%   [im, C_co] = COSEM_OSL(im, D, beta, dU, A, uu, epps, C_co, h, COSEM_MAP, osa_iter, SinDelayed, is_transposed)
%   [im, C_co] = COSEM_OSL(im, D, beta, dU, A, uu, epps, C_co, h, COSEM_MAP, osa_iter, SinDelayed, is_transposed, options, Nx, Ny, Nz, gaussK)
%   [im, C_co] = COSEM_OSL(im, D, beta, dU, RHS, osa_iter, h, C_co, COSEM_MAP)
% INPUTS:
%   im = The current estimate
%   D = Sum of the complete data system matrix (D = sum(B,2), where B =
%   [A_1,A_2,...,A_subsets]
%   beta = Regularization parameter
%   dU = Gradient of the prior
%   A = The transpose of the (sparse) system matrix at current subset
%   uu = Measurements at current subset
%   epps = Small constant to prevent division by zero
%   C_co = Complete-data matrix
%   h = Acceleration parameter (ACOSEM)
%   COSEM_MAP = Whether COSEM or ACOSEM is used (1 == ACOSEM)
%   osa_iter = Current subset (sub-iteration)
%   SinDelayed = Randoms and/or scatter correction data. Dimension must be
%   either a scalar or a vector of same size as uu. If no scatter and/or
%   randoms data is available, use zero.
%   is_transposed = true if A matrix is the transpose of it, false if not
%
%   These values are needed only when doing PSF reconstruction:
%	options = Needed for PSF reconstruction
%   Nx = Image size in X-direction, for PSF
%   Ny = Image size in Y-direction, for PSF
%   Nz = Image size in Z-direction, for PSF
%   gaussK = Gaussian kernel for PSF
%
%   This is for implementation 4:
%   RHS = Backprojection computed by implementation 4
%
% OUTPUTS:
%   im = The updated estimate
%   C_co = Updated Complete-data matrix
%
%   See also COSEM_im, ACOSEM_im, MBSREM, OSL_OSEM, BSREM_iter, RBI_subiter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ville-Veikko Wettenhovi
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
% Implementation 4
if isempty(osa_iter) && isempty(is_transposed)
    if h == 1
        C_co(:,uu) = (im).^(1 / epps) .* A;
        im = (sum(C_co,2) ./ (D + beta * dU)) .^epps;
    else
        C_co(:,uu) = (im) .* A;
        im = (sum(C_co,2) ./ (D + beta * dU));
    end
else
    % Implementation 1
    % PSF reconstruction
    if ~isempty(varargin{1}) && varargin{1}.use_psf
        im_apu = computeConvolution(im, varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5});
    else
        im_apu = im;
    end
    if is_transposed
        FP = A' * im_apu + epps + SinDelayed;
        if ~isempty(varargin{1}) && isfield(varargin{1},'CT') && varargin{1}.CT
            FP = exp(FP);
        end
        RHS = A * (uu ./ FP);
    else
        FP = A * im_apu + epps + SinDelayed;
        if ~isempty(varargin{1}) && isfield(varargin{1},'CT') && varargin{1}.CT
            FP = exp(FP);
        end
        RHS = A' * (uu ./ FP);
    end
    if ~isempty(varargin) && ~isempty(varargin{1}) && varargin{1}.use_psf
        RHS = computeConvolution(RHS, varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5});
    end
    % ACOSEM-OSL
    if COSEM_MAP == 1
        C_co(:,osa_iter) = im.^(1/h) .* RHS;
        im = (sum(C_co,2)./(D + beta * dU)).^h;
        if ~isempty(varargin{1}) && varargin{1}.use_psf
            im_apu = computeConvolution(im, varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5});
        else
            im_apu = im;
        end
        if is_transposed
            if ~isempty(varargin) && ~isempty(varargin{1}) && isfield(varargin{1}, 'CT') && varargin{1}.CT
                im = (im)*(sum(uu)/(sum(exp(A'*im_apu))));
            else
                im = (im)*(sum(uu)/(sum(A'*im_apu)));
            end
        else
            if ~isempty(varargin) && ~isempty(varargin{1}) && isfield(varargin{1}, 'CT') && varargin{1}.CT
                im = (im)*(sum(uu)/(sum(exp(A*im_apu))));
            else
                im = (im)*(sum(uu)/(sum(A*im_apu)));
            end
        end
    else
        % COSEM-OSL
        C_co(:,osa_iter) = im .* RHS;
        im = (sum(C_co,2)./(D + beta * dU));
    end
    im(im < epps) = epps;
end