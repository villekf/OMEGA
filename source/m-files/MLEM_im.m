function [im, varargout] = MLEM_im(im, Summ, epps, varargin)
%MLEM_im Computes the Maximum Likelihood Expectation Maximization (MLEM)
%estimates
%
% Example:
%   im = MLEM_im(im, Summ, epps, A, uu, is_transposed)
% INPUTS:
%   im = The current estimate
%   Summ = sum(A,2) when is_transposed = true and sum(A,1)' when false
%   (sensitivity image)
%   epps = Small constant to prevent division by zero
%   A = The (transpose of the) (sparse) system matrix at current subset
%   uu = Measurement data
%   is_transposed = true if A matrix is the transpose of it, false if not
%   options = use_psf is used
%   Nx = Image size in x-direction
%   Ny = y-direction
%   Nz = z-direction
%   gaussK = PSF kernel
%
% OUTPUTS:
%   im = The updated estimate
%   log_lh = The log-likelihood value
% 
% See also OSEM_im, ROSEM_subiter, BSREM_subiter, DRAMA_subiter,
% RBI_subiter, COSEM_im, ACOSEM_im, ECOSEM_im, MBSREM

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
if nargin >= 6
    if ~isempty(varargin) && length(varargin) > 3 && ~isempty(varargin{4}) && varargin{4}.use_psf
        im_apu = computeConvolution(im, varargin{4}, varargin{5}, varargin{6}, varargin{7}, varargin{8});
    else
        im_apu = im;
    end
    if varargin{3}
        if length(varargin) > 3 && ~isempty(varargin{4}) && isfield(varargin{4},'CT') && varargin{4}.CT
            RHS = (varargin{1} * (varargin{2} ./ (exp(varargin{1}' * im_apu))));
        else
            RHS = (varargin{1} * (varargin{2} ./ (varargin{1}' * im_apu + epps)));
        end
    else
        if length(varargin) > 3 && ~isempty(varargin{4}) && isfield(varargin{4},'CT') && varargin{4}.CT
            RHS = (varargin{1}' * (varargin{2} ./ (exp(varargin{1} * im_apu))));
        else
            RHS = (varargin{1}' * (varargin{2} ./ (varargin{1} * im_apu + epps)));
        end
    end
    if ~isempty(varargin) && length(varargin) > 3 && ~isempty(varargin{4}) && varargin{4}.use_psf
        RHS = computeConvolution(im, varargin{4}, varargin{5}, varargin{6}, varargin{7}, varargin{8});
    end
    im = (im ./ Summ) .* RHS;
    if nargout > 1
        if varargin{3}
            varargout{1} = sum(varargin{2}.*log(varargin{1}' * im)-(varargin{1}' * im)-log(factorial(round(varargin{2}))));
        else
            varargout{1} = sum(varargin{2}.*log(varargin{1} * im)-(varargin{1} * im)-log(factorial(round(varargin{2}))));
        end
    end
elseif nargin == 4
    im = (im ./ Summ) .* (varargin{1} + epps);
else
    error('Invalid number of input arguments')
end