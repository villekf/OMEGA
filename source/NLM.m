function output = NLM(input,Ndx, Ndy, Ndz, Nlx, Nly, Nlz,h2, epps, Nx, Ny, Nz, options)
% Non-Local Means prior (NLM)
% Based on Simple Non Local Means (NLM) Filter from MATLAB file exchange 
% https://se.mathworks.com/matlabcentral/fileexchange/52018-simple-non-local-means-nlm-filter
% 
% INPUTS:
% im = The current estimate
% Ndx = Distance of search window in X-direction
% Ndy = Distance of search window in Y-direction
% Ndz = Distance of search window in Z-direction
% Nlx = Patch size in X-direction
% Nly = Patch size in Y-direction
% Nlz = Patch size in Z-direction
% h2 = NLM filter parameter
% epps = Small constant to prevent division by zero
% Nx = Image (estimate) size in X-direction
% Ny = Image (estimate) size in Y-direction
% Nz = Image (estimate) size in Z-direction
% options = Reconstruction options (anatomical, MRP algorithm)
%
% OUTPUTS:
% output = The (gradient of) NLM prior
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi
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

input = reshape(input, Nx, Ny, Nz);
[m, n, k]=size(input);
pixels = input(:);

s = m*n*k;

psizex = 2*Nlx+1;
psizey = 2*Nly+1;
psizez = 2*Nlz+1;
nsizex = 2*Ndx+1;
nsizey = 2*Ndy+1;
nsizez = 2*Ndz+1;


% Compute list of edges (pixel pairs within the same search window)
indexes = reshape(1:s, m, n, k);
padIndexes = padding(indexes, [Ndx Ndy Ndz]);
neighbors = im2col_3D_sliding(padIndexes, [nsizex nsizey nsizez]);
neighbors = reshape(neighbors, [], size(neighbors,2)*size(neighbors,3));
TT = repmat(1:s, [nsizex*nsizey*nsizez 1]);
edges = [TT(:) neighbors(:)];
RR = TT(:) >= neighbors(:);
edges(RR, :) = [];

% Compute weight matrix (using weighted Euclidean distance)
if options.NLM_use_anatomical
    padInput = padding(options.NLM_ref,[Nlx Nly Nlz]);
    patches = im2col_3D_sliding(padInput, [psizex psizey psizez]);
    patches = reshape(patches, [], size(patches,2)*size(patches,3))';
    diff = patches(edges(:,1), :) - patches(edges(:,2), :);
else
    % Compute patches
    %  padInput = padarray(input,[Nlx Nly Nlz],'symmetric');
    padInput = padding(input,[Nlx Nly Nlz]);
    patches = im2col_3D_sliding(padInput, [psizex psizey psizez]);
    patches = reshape(patches, [], size(patches,2)*size(patches,3))';
    diff = patches(edges(:,1), :) - patches(edges(:,2), :);
end
V = exp(-sum(diff.*diff,2)/h2^2);
V(V == 0) = epps;
W = sparse(edges(:,1), edges(:,2), V, s, s);


W = W + W';

% Normalize weights
W = spdiags(1./(sum(W,2)), 0, s, s)*W;

if options.NLM_MRP
    output = pixels - W*pixels;
else
    S = pixels(edges(:,1), :) - pixels(edges(:,2), :);
    
    K = sparse(edges(:,1), edges(:,2), S, s, s);
    K = K + K';
    
    
    output = full(sum(W .* K,2));
end

% NLTV
%  output = full(sum(W .* K,2) ./ (sqrt(sum(W .* K.^2,2)) + epps));