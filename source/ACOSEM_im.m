function [im, C_aco] = ACOSEM_im(im, A, epps, uu, C_aco, D, h, osa_iter, SinDelayed, is_transposed, varargin)
%ACOSEM_IM Computes the Accelerated COSEM (ACOSEM) estimates
%
% Examples:
%   [im, C_aco] = ACOSEM_im(im, A, epps, uu, C_aco, D, h, osa_iter, SinDelayed, is_transposed)
%   [im, C_aco] = ACOSEM_im(im, A, epps, uu, C_aco, D, h, osa_iter, SinDelayed, is_transposed, options, Nx, Ny, Nz, gaussK)
%   [im, C_aco] = ACOSEM_im(im, RHS, C_aco, D, osa_iter, h)
% INPUTS:
%   im = The current estimate
%   A = The transpose of the (sparse) system matrix at current subset
%   epps = Small constant to prevent division by zero
%   uu = Measurements at current subset
%   C_aco = Complete-data matrix for ACOSEM
%   D = Sum of the complete data system matrix (D = sum(B,2), where B =
%   [A_1,A_2,...,A_subsets]
%   h = Acceleration parameter
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
%   C_aco = Updated complete-data matrix
%
%   See also COSEM_im, ECOSEM_im, COSEM_OSL

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
% This section is for implementation 4
if isempty(osa_iter) && isempty(is_transposed)
    % Compute the scaled ACOSEM estimate
    if isempty(uu)
        im = (im)*(sum(epps)/(sum(A)));
        C_aco = [];
    else
        % Compute the non-scaled ACOSEM estimate
        epps(:,C_aco) = im.^(1/D) .* A;
        im = (sum(epps,2) ./ uu).^D;
        C_aco = epps;
    end
else
    % PSF blurring
    if ~isempty(varargin) && ~isempty(varargin{1}) && varargin{1}.use_psf
        im_apu = computeConvolution(im, varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5});
    else
        im_apu = im;
    end
    if is_transposed
        FP = A' * im_apu + epps + SinDelayed;
        % Backprojection
        RHS = A * (uu ./ FP);
    else
        FP = A * im_apu + epps + SinDelayed;
        RHS = A' * (uu ./ FP);
    end
    if ~isempty(varargin) && ~isempty(varargin{1}) && varargin{1}.use_psf
        RHS = computeConvolution(RHS, varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5});
    end
    C_aco(:,osa_iter) = im.^(1/h) .* RHS;
    im = (sum(C_aco,2)./D).^h;
    if ~isempty(varargin) && ~isempty(varargin{1}) && varargin{1}.use_psf
        im_apu = computeConvolution(im, varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5});
    else
        im_apu = im;
    end
    if is_transposed
        im = (im)*(sum(uu)/(sum(A'*im_apu)));
    else
        im = (im)*(sum(uu)/(sum(A*im_apu)));
    end
    im(im < epps) = epps;
end