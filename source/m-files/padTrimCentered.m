function B = padTrimCentered(A, n, dim)
%PADTRIMCENTERED Pad with zeros or trim an array symmetrically along one dimension
%   B = PADTRIMCENTERED(A, N, DIM) pads with zeros or trims the array A
%   along dimension DIM so that size(B, DIM) == N. With d = abs(N -
%   size(A, DIM)), the leading edge gets floor(d/2) and the trailing edge
%   ceil(d/2) elements, both when padding and when trimming. The class of
%   A is preserved. This reproduces MATLAB resize(A, N, Dimension=DIM,
%   FillValue=0, Side="both") without requiring that function.
%
% See also projectorClass

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2026 Ville-Veikko Wettenhovi
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

sz = size(A);
if dim > numel(sz)
    sz(end + 1 : dim) = 1;
end
m = sz(dim);
if n == m
    B = A;
    return;
end
d = abs(n - m);
lead = floor(d / 2);
subs = repmat({':'}, 1, numel(sz));
if n > m
    szB = sz;
    szB(dim) = n;
    B = zeros(szB, class(A));
    subs{dim} = lead + 1 : lead + m;
    B(subs{:}) = A;
else
    subs{dim} = lead + 1 : lead + n;
    B = A(subs{:});
end
end
