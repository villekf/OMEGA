function s = powerMethod(A, varargin)
%POWERMETHOD Computes the power method for the largest singular value in a
%matrix
%   Examples:
%       s = powerMethod(A)
%       s = powerMethod(A, nIter)
%       s = powerMethod(A, nIter, useTV)
%       s = powerMethod(A, nIter, useTV, Nx, Ny, Nz)
%   Inputs:
%       A = Either the system matrix or the forward/backward projection
%       class object.
%       nIter = (Optional) Number of iterations, default value is 10.
%       useTV = (Optional) Whether a TV CP-algorithm is used. Default is
%       false (no TV).
%       Nx = The number of voxels in the x-direction (required only when A
%       is a matrix)
%       Ny = The number of voxels in the y-direction (required only when A
%       is a matrix)
%       Nz = The number of voxels in the z-direction (required only when A
%       is a matrix)
%   Output:
%       s = The estimate for the largest singular value.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021-2024 Ville-Veikko Wettenhovi
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
if nargin >= 2 || isempty(varargin) || isempty(varargin{1})
    Niter = 10;
else
    Niter = varargin{1};
end
% if nargin >= 3 || isempty(varargin) || isempty(varargin{2})
%     useTV = false;
% else
%     useTV = varargin{2};
% end
if nargin >= 3|| isempty(varargin) || isempty(varargin{2})
    useFilt = false;
else
    useFilt = true;
    Nf = 2^nextpow2(varargin{2}(1));
    filt = rampFilt(Nf);
end
if isa(A,'projectorClass')
    useClass = true;
else
    useClass = false;
end
% x0 = ones(dim,'single')*0.00002;
if useClass
    x = cell(A.param.nMultiVolumes + 1, 1);
    nVol = A.param.nMultiVolumes + 1;
    for kk = 0 : A.param.nMultiVolumes
        if (kk > 0 && mod(kk,2) == 0)
            x{kk + 1} = x{kk};
        else
            x{kk + 1} = abs(randn(A.param.N(kk + 1), 1, A.param.cType));
            x{kk + 1} = x{kk + 1} ./ norm(x{kk + 1},2);
        end
    end
    Niter = A.param.powerIterations;
    s = zeros(nVol, 1, A.param.cType);
    A.param.subset = 1;
else
    dim = size(A,2);
    x = randn(dim, 1);
    x = x ./ norm(x,2);
    nVol = 1;
end
if useClass
    if A.param.implementation == 1
        if A.param.subsets > 1
            B = formMatrix(A, 1, 1);
        else
            B = formMatrix(A, 0, 1);
        end
    end
end
for kk = 1 : Niter
    if useClass && A.param.nMultiVolumes > 0
        if A.param.implementation == 1
            temp = B' * x;
        else
            temp = A * x;
        end
    else
        if useClass && A.param.implementation == 1
            temp = B' * x{1};
        else
            temp = A * x{1};
        end
    end
    if useFilt
        temp = reshape(temp, varargin{6});
        padf = zeros(Nf, size(temp,2), size(temp,3));
        padf(Nf/2 - size(temp,1)/2+1:Nf/2 +size(temp,1)/2,:,:) = temp;
        ft = real(ifft(fft(padf) .* filt));
        temp = ft(Nf/2 - size(temp,1)/2+1:Nf/2 +size(temp,1)/2,:,:);
        temp = temp(:);
    elseif useClass
        temp = applyMeasPreconditioning(A.param, temp, A.nMeasSubset(1), 1);
    end
    if useClass && A.param.implementation == 1
        temp = B * temp;
    else
        temp = A' * temp;
    end
    for ii = 1 : nVol
        if useClass
            if A.param.nMultiVolumes > 0
                temp{ii} = applyImagePreconditioning(A.param, temp{ii}, x{ii}, kk, ii);
                s(ii) = (x{ii}' * temp{ii} * A.param.subsets) ./ (x{ii}' * x{ii});
                x{ii} = temp{ii};
            else
                temp = applyImagePreconditioning(A.param, temp, x{ii}, kk, ii);
                s(ii) = (x{ii}' * temp * A.param.subsets) ./ (x{ii}' * x{ii});
                x{ii} = temp;
            end
            x{ii} = x{ii} ./ norm(x{ii});
            if A.param.verbose >= 2
                disp(['Largest eigenvalue at iteration ' num2str(kk) ' is ' num2str(s(ii))])
            end
        else
            s = (x' * temp) ./ (x' * x);
            x = temp;
            x = x ./ norm(x);
        end
    end
end
