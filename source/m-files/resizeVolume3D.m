function B = resizeVolume3D(A, sz)
%RESIZEVOLUME3D Resize a 3-D volume to the size SZ
%   B = RESIZEVOLUME3D(A, SZ) resizes the volume A to the size SZ = [n1 n2
%   n3]. If imresize3 is available it is used (default interpolation). If
%   not (e.g. GNU Octave), a trilinear fallback with the same pixel-centre
%   convention as imresize is used: the output pixel centres map to the
%   input coordinates (c_out - 0.5) * size(A,d) / SZ(d) + 0.5, clamped to
%   the input volume. The fallback returns single if A is single and double
%   otherwise.
%
% See also imresize3, interpn

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
if exist('imresize3', 'file') == 2
    B = imresize3(A, sz);
    return
end

sz = sz(:)';
if numel(sz) < 3
    sz(end + 1 : 3) = 1;
end
isSingle = isa(A, 'single');
A = double(A);
c = cell(1, 3);
for d = 1 : 3
    nIn = size(A, d);
    ci = ((1 : sz(d)) - 0.5) * nIn / sz(d) + 0.5;
    c{d} = min(max(ci, 1), nIn);
    if nIn == 1
        % interpn needs at least two points per dimension
        A = cat(d, A, A);
    end
end
[C1, C2, C3] = ndgrid(c{1}, c{2}, c{3});
B = interpn(A, C1, C2, C3, 'linear');
if isSingle
    B = single(B);
end
end
