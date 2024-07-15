function [grad, varargout] = Weighted_mean(im, weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, mean_type, varargin)
%Weighted_mean Weighted mean prior
%   This function computes the gradient of the weighted mean prior. Two
%   different prior types are available that each can compute three
%   different means (arithmetic, harmonic, geometric). The first three mean
%   types (mean_type = 1,2,3) compute the means in an MRP fashion, e.g. the
%   median is simply replaced by the weighted mean. The last three
%   (mean_type = 4,5,6) compute a prior that computes the mean around the
%   neighborhood of each pixel. The latter three also have differentiable
%   MRF-functions and allows the output of the MRF-functions with the
%   specified input parameters (second output value).
%
% Examples:
%   grad = Weighted_mean(im, weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, mean_type)
%   grad = Weighted_mean(im, weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, mean_type, epps, med_no_norm)
% INPUTS:
%   im = The current estimate
%   weighted_weights = Weights for the weighted mean
%   Nx = Image (estimate) size in X-direction
%   Ny = Image (estimate) size in Y-direction
%   Nz = Image (estimate) size in Z-direction
%   Ndx = How many pixels included in X-direction neighborhood
%   Ndy = How many pixels included in Y-direction neighborhood
%   Ndz = How many pixels included in Z-direction neighborhood
%   mean_type = Type of the mean, 1 = Arithmetic MRP, 2 = Harmonic MRP, 3 =
%   Geometric MRP, 4 = Arithmetic prior, 5 = Harmonic prior, 6 = Geometric
%   prior
%   epps = Small constant to prevent division by zero, needed only if
%   med_no_norm = true (can be omitted)
%   med_no_norm = if true, then no normalization will be performed in the
%   MRP-based priors (division by the weighted mean filtered image) (can be
%   omitted, assumed to be false then)
%
% OUTPUTS:
%   grad = The (gradient of) weighted mean prior
%
% See also MRP, TGV, L_filter, Quadratic_prior, TVpriorfinal, FMH

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
if Nz==1
    padd = padding(reshape(im,Nx,Ny),[Ndx Ndy]);
else
    padd = padding(reshape(im,Nx,Ny,Nz),[Ndx Ndy Ndz]);
end
if nargin >= 10
    if ~isempty(varargin{1})
        epps = varargin{1};
    else
        epps = 1e-8;
    end
    med_no_norm = varargin{2};
else
    med_no_norm = false;
end
w_sum = sum(weighted_weights(:));
if mean_type == 1
    grad = convn(padd, weighted_weights ./ w_sum, 'valid');
elseif mean_type == 2
    grad = 1 ./ convn(1 ./ padd, weighted_weights ./ w_sum, 'valid');
elseif mean_type == 3
    grad = exp(convn(log(padd), weighted_weights ./ w_sum, 'valid'));
elseif mean_type == 4
    if Nz==1
        padd = padding(reshape(im,Nx,Ny),[Ndx*2 Ndy*2]);
    else
        padd = padding(reshape(im,Nx,Ny,Nz),[Ndx*2 Ndy*2 Ndz*2]);
    end
    m = convn(padd, weighted_weights ./ w_sum, 'valid');
    mm = zeros(size(im,1),numel(weighted_weights));
    jj = 1;
    for ll = 1 : size(weighted_weights,3)
        for kk = 1: size(weighted_weights,2)
            for hh = 1 : size(weighted_weights,1)
                apu = m(hh: hh + Nx - 1, kk: kk + Ny - 1, ll : ll + Nz - 1);
                mm(:,jj) = apu(:);
                jj = jj + 1;
            end
        end
    end
    grad = bsxfun(@rdivide, im, mm);
    grad(grad < 1e-3) = 1e-3;
    grad = bsxfun(@times, log(grad), weighted_weights(:)');
    grad = sum(grad,2);
    grad = grad(:);
    if nargout >= 2
        mrf = bsxfun(@minus, bsxfun(@times, im, log(bsxfun(@rdivide, im, mm))), im) + mm;
        mrf = bsxfun(@times, mrf, weighted_weights(:)');
        mrf = sum(mrf,2);
        varargout{1} = mrf;
    end
elseif mean_type == 5
    if Nz==1
        padd = padding(reshape(im,Nx,Ny),[Ndx*2 Ndy*2]);
    else
        padd = padding(reshape(im,Nx,Ny,Nz),[Ndx*2 Ndy*2 Ndz*2]);
    end
    m = 1./convn(1 ./ padd, weighted_weights ./ w_sum, 'valid');
    mm = zeros(size(im,1),numel(weighted_weights));
    jj = 1;
    for ll = 1 : size(weighted_weights,3)
        for kk = 1: size(weighted_weights,2)
            for hh = 1 : size(weighted_weights,1)
                apu = m(hh: hh + Nx - 1, kk: kk + Ny - 1, ll : ll + Nz - 1);
                mm(:,jj) = apu(:);
                jj = jj + 1;
            end
        end
    end
    im(im < 1e-1) = 1e-1;
    grad = bsxfun(@rdivide, bsxfun(@minus, im, mm), im) - bsxfun(@rdivide, bsxfun(@minus, im, mm).^2, 2*im.^2);
    grad = bsxfun(@times, grad, weighted_weights(:)');
    grad = sum(grad,2);
    grad = grad(:);
    if nargout >= 2
        mrf = bsxfun(@rdivide, bsxfun(@minus, im, f), 2*im);
        mrf = bsxfun(@times, mrf, weighted_weights(:)');
        mrf = sum(mrf,2);
        varargout{1} = mrf;
    end
elseif mean_type == 6
    if Nz==1
        padd = padding(reshape(im,Nx,Ny),[Ndx*2 Ndy*2]);
    else
        padd = padding(reshape(im,Nx,Ny,Nz),[Ndx*2 Ndy*2 Ndz*2]);
    end
    m = exp(convn(log(padd), weighted_weights ./ w_sum, 'valid'));
    mm = zeros(size(im,1),numel(weighted_weights));
    jj = 1;
    for ll = 1 : size(weighted_weights,3)
        for kk = 1: size(weighted_weights,2)
            for hh = 1 : size(weighted_weights,1)
                apu = m(hh: hh + Nx - 1, kk: kk + Ny - 1, ll : ll + Nz - 1);
                mm(:,jj) = apu(:);
                jj = jj + 1;
            end
        end
    end
    im(im < 1e-1) = 1e-1;
    grad = 1 - bsxfun(@rdivide, mm, im);
    grad = bsxfun(@times, grad, weighted_weights(:)');
    grad = sum(grad,2);
    grad = grad(:);
    if nargout >= 2
        mrf = bsxfun(@plus, bsxfun(@times, mm, log(bsxfun(@rdivide, mm, im))), im) - mm;
        mrf = bsxfun(@times, mrf, weighted_weights(:)');
        mrf = sum(mrf,2);
        varargout{1} = mrf;
    end
else
    error('Unsupported mean')
end
if mean_type <= 3
    if med_no_norm
        grad = (im - grad(:));
    else
        grad = (im - grad(:))./(grad(:) + epps);
    end
end