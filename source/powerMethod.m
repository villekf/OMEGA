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
if isempty(varargin) || isempty(varargin{1})
    Niter = 10;
else
    Niter = varargin{1};
end
if isempty(varargin) || isempty(varargin{2})
    useTV = false;
else
    useTV = varargin{2};
end
if nargin >= 6 && (~isempty(varargin) || (~isempty(varargin{3}) && ~isempty(varargin{4}) && ~isempty(varargin{5})))
    dim = [varargin{3} varargin{4} varargin{5}];
else
    if ismatrix(A) && ~isa(A,'forwardBackwardProjectCT') && ~isa(A,'forwardBackwardProject')
        error('When A is a matrix, the image size needs to be input')
    else
        dim = [A.OProperties.Nx A.OProperties.Ny A.OProperties.Nz];
    end
end
x0 = ones(dim)*0.02;
x = x0(:);
N = prod(dim);
if useTV
    gradx = zeros(dim);
    grady = zeros(dim);
    gradz = zeros(dim);
    divx = zeros(dim);
    divy = zeros(dim);
    divz = zeros(dim);
end
for kk = 1 : Niter
    temp = A * x;
    temp = A' * temp;
    if useTV
        gr = reshape(x,dim);
        gradx(:,1:dim(2)-1,:) = diff(gr,1,2);
        grady(1:dim(1)-1,:,:) = diff(gr,1,1);
        if dim(3) > 1
            gradz(:,:,1:dim(3)-1) = diff(gr,1,3);
            grad = [gradx(:);grady(:);gradz(:)];
        else
            grad = [gradx(:);grady(:)];
        end
        qx = reshape(grad(1:N),dim);
        qy = reshape(grad(N+1:2*N),dim);
        if dim(3) > 1
            qz = reshape(grad(N*2+1:end),dim);
        end
        divx(:,2:dim(2),:) = -diff(qx,1,2);
        divx(:,1,:) = -qx(:,1,:);
        divx(:,dim(2),:) = qx(:,dim(2)-1,:);
        divy(2:dim(1),:,:) = -diff(qy,1,1);
        divy(1,:,:) = -qy(1,:,:);
        divy(dim(1),:,:) = qy(dim(1)-1,:,:);
        if dim(3) > 1
            divy(:,:,2:dim(3)) = -diff(qz,1,3);
            divy(:,:,1) = -qz(:,:,1);
            divy(:,:,dim(3)) = qz(:,:,dim(3)-1);
            grad = (divx+divy+divz);
        else
            grad = (divx+divy);
        end
        x = temp + grad(:);
        x = x ./ norm(x,2);
    else
        x = temp;
        x = x ./ norm(x,2);
    end
end
if useTV
    s = sqrt(norm(A*x,2).^2 + norm(grad(:),2).^2);
else
    s = norm(A*x,2);
end

