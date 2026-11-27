function B = resampleDetectorAxis(A, dim, nOut, ratio)
%RESAMPLEDETECTORAXIS Resamples one axis of a detector image onto a lattice of another pixel pitch, centered on the detector
%   B = resampleDetectorAxis(A, dim, nOut, ratio) resamples the dimension dim (1 or 2) of A to nOut pixels. ratio is the pitch of the
%   output pixels divided by the pitch of the input pixels (not necessarily an integer), the centers of both lattices are at the
%   center of the input array. The output pixel j is at the input index  u(j) = (j - (nOut + 1) / 2) * ratio + (size(A, dim) + 1) / 2.
%
%   The interpolation is that of imresize with its default (bicubic, i.e. Keys cubic convolution, a = -0.5) kernel, including the
%   antialiasing (kernel stretched by the ratio) when ratio > 1: the output is a weighted mean of the input values (the weights of
%   a pixel sum to one), so the values and not the sums of the counts are preserved, exactly like imresize. Input samples outside
%   the array are zeros (the weights are normalized before they are dropped, so the edge of the array is dark, like with zero padding
%   followed by imresize). For an integer ratio of a symmetric window this is the same as padding or trimming with resize followed by
%   imresize, except near the edges of a trimmed window (where imresize mirrors the window, here the real neighbouring pixels are used).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021-2025 Ville-Veikko Wettenhovi
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nIn = size(A, dim);
nOut = double(nOut);
ratio = double(ratio);

if ratio > 1 % Shrinking: antialiasing, as in imresize
    h = @(x) (1 / ratio) * cubicKernel(x / ratio);
    kernelWidth = 4 * ratio;
else
    h = @cubicKernel;
    kernelWidth = 4;
end
j = (1 : nOut)';
u = (j - (nOut + 1) / 2) * ratio + (nIn + 1) / 2;
left = floor(u - kernelWidth / 2);
P = ceil(kernelWidth) + 2;
indices = left + (0 : P - 1);
weights = h(u - indices);
weights = weights ./ sum(weights, 2);
weights(indices < 1 | indices > nIn) = 0; % Zero outside of the detector
W = zeros(nOut, nIn);
for kk = 1 : P
    ok = indices(:, kk) >= 1 & indices(:, kk) <= nIn;
    W(sub2ind(size(W), find(ok), indices(ok, kk))) = W(sub2ind(size(W), find(ok), indices(ok, kk))) + weights(ok, kk);
end

cls = class(A);
if ~isfloat(A)
    A = double(A);
end
W = cast(W, 'like', A);
sz = size(A);
sz(end + 1 : max(dim, 2)) = 1;
if dim == 1
    B = reshape(W * reshape(A, sz(1), []), [nOut, sz(2 : end)]);
else
    perm = [2, 1, 3 : numel(sz)];
    B = permute(reshape(W * reshape(permute(A, perm), sz(2), []), [nOut, sz(perm(2 : end))]), perm);
end
if ~strcmp(cls, class(B))
    B = cast(B, cls);
end
end

function y = cubicKernel(x)
% Cubic convolution kernel (a = -0.5), the default of imresize
absx = abs(x);
absx2 = absx .^ 2;
absx3 = absx .^ 3;
y = (1.5 * absx3 - 2.5 * absx2 + 1) .* (absx <= 1) + (-0.5 * absx3 + 2.5 * absx2 - 4 * absx + 2) .* ((1 < absx) & (absx <= 2));
end
