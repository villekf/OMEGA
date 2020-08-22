function output = NLM(input,Ndx, Ndy, Ndz, Nlx, Nly, Nlz, h2, epps, Nx, Ny, Nz, options)
% Non-Local Means prior (NLM)
%
%   Supports regular NLM, NLTV and NLM-MRP regularization. Also allows the
%   weight matrix to be determined from a reference image.
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
% options = Reconstruction options (anatomical, NLTV, MRP algorithm)
%
% OUTPUTS:
% output = The (gradient of) NLM prior
%

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

input = reshape(input, Nx, Ny, Nz);
padx = Nlx + Ndx;
pady = Nly + Ndy;
padz = Nlz + Ndz;
% Padd the images (symmetric)
if options.NLM_use_anatomical
    padInput = padding(options.NLM_ref,[padx pady padz]);
else
    padInput = padding(input,[padx pady padz]);
end
padInput = padInput(:);
input = padding(input,[padx pady padz]);
[N, M, K] = size(input);
N = uint32(N);
M = uint32(M);
K = uint32(K);
input = input(:);
% What type of NLM
if options.NLTV
    type = int32(1);
elseif options.NLM_MRP
    type = int32(2);
else
    type = int32(0);
end
if exist('OCTAVE_VERSION','builtin') == 0
%%% MEX
    output = NLM_func(padInput, input, options.gaussianNLM, int32(Ndx), int32(Ndy), int32(Ndz), int32(Nlx), int32(Nly), int32(Nlz),...
        N, M, K, h2*h2, type, epps);
elseif exist('OCTAVE_VERSION','builtin') == 5
%%% OCT
    output = NLM_oct(padInput, input, options.gaussianNLM, int32(Ndx), int32(Ndy), int32(Ndz), int32(Nlx), int32(Nly), int32(Nlz),...
        N, M, K, h2*h2, type, epps);
end
output = reshape(output, N, M, K);
% Convert back to original image size
output = output(padx + 1 : end - padx,pady + 1 : end - padx,padz + 1 : end - padz);
output = output(:);
if options.NLM_MRP && type == 2
    input = input(padx + 1 : end - padx,pady + 1 : end - padx,padz + 1 : end - padz);
    output = input - output;
end
end
